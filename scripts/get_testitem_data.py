"""Command line interface for retrieving ChEMBL and PubChem compound data."""

from __future__ import annotations

import json
import sys

# ruff: noqa: E402
from pathlib import Path

if __package__ is None:  # running as a script
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import argparse
from collections import ChainMap
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from itertools import islice

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library import molecule_catalog
from library import testitem_enrichment
from library import pubchem_library as pl
from library.molecule_catalog import (
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    ApiCfg,
    Config,
    MoleculeCatalogCfg,
    PubChemCfg,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_testitems
from schemas import TestitemsSchema, normalize_testitems

_TYPO_PARENT_COLUMN = "parant_molecule_id"


def ensure_no_parant_column(df: pd.DataFrame) -> None:
    """Raise a :class:`ValueError` if the legacy typo column is present."""

    if _TYPO_PARENT_COLUMN in df.columns:
        raise ValueError(
            "unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead"
        )


PARENT_LOOKUP_SOURCE_CACHE = "cache"
PARENT_LOOKUP_SOURCE_REMOTE = "remote"
PARENT_LOOKUP_SOURCE_SKIPPED = "skipped"


PUBCHEM_COLUMNS = [
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
]

PUBCHEM_CID_CACHE_FILENAME = "pubchem_cid_cache.json"
_CACHE_MISSING = object()


@dataclass(frozen=True)
class ParentLookupStats:
    """Summary information about parent molecule enrichment."""

    source: str
    missing: int
    unique: int
    attached: int


def _cache_state(path: Path) -> tuple[bool, float | None]:
    """Return cache file presence and modification time."""

    try:
        stat = path.stat()
    except FileNotFoundError:
        return False, None
    return True, stat.st_mtime


def _resolve_parent_source(
    before: tuple[bool, float | None], after: tuple[bool, float | None]
) -> str:
    """Determine parent lookup source from cache state transitions."""

    before_exists, before_mtime = before
    after_exists, after_mtime = after

    if after_exists and before_exists and after_mtime == before_mtime:
        return PARENT_LOOKUP_SOURCE_CACHE
    if after_exists or before_exists != after_exists:
        return PARENT_LOOKUP_SOURCE_REMOTE
    return PARENT_LOOKUP_SOURCE_REMOTE


def _normalise_chembl_ids(series: pd.Series) -> pd.Series:
    """Return ``series`` normalised to upper-case ChEMBL identifiers."""

    normalised = (
        series.fillna("")
        .astype("string")
        .str.strip()
        .str.upper()
    )
    return normalised


def _normalise_identifier(value: object, *, uppercase: bool = False) -> str | None:
    """Return a stripped string representation or ``None`` for missing values."""

    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    if isinstance(value, pd.Series):  # pragma: no cover - defensive
        raise TypeError("_normalise_identifier expects scalar values")
    if pd.isna(value):  # handles pd.NA, numpy.nan
        return None
    text = str(value).strip()
    if not text:
        return None
    return text.upper() if uppercase else text


def _load_pubchem_cid_cache(path: Path) -> dict[str, str | None]:
    """Load a CID cache from ``path`` if it exists."""

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return {}
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("pubchem_cache_load_failed", path=str(path), error=str(exc))
        return {}

    cache: dict[str, str | None] = {}
    for key, raw_value in data.items():
        chembl_id = _normalise_identifier(key, uppercase=True)
        if not chembl_id:
            continue
        cid = (
            _normalise_identifier(raw_value)
            if isinstance(raw_value, str)
            else None
        )
        cache[chembl_id] = cid
    return cache


def _dump_pubchem_cid_cache(cache: Mapping[str, str | None], path: Path) -> None:
    """Persist ``cache`` to ``path`` in JSON format."""

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        serialisable = {
            key: value for key, value in sorted(cache.items()) if key
        }
        path.write_text(
            json.dumps(serialisable, indent=2, sort_keys=True),
            encoding="utf-8",
        )
    except OSError as exc:
        logger.warning("pubchem_cache_write_failed", path=str(path), error=str(exc))


def _select_primary_cid(
    cid_response: str | None,
    *,
    chembl_id: str | None,
    identifier: str,
    source: str,
) -> str | None:
    """Return the first CID from ``cid_response`` and log additional matches."""

    if not cid_response:
        return None
    candidates = [c for c in cid_response.split("|") if c]
    if not candidates:
        return None
    primary = candidates[0]
    if len(candidates) > 1:
        logger.info(
            "pubchem_multiple_cids",
            chembl_id=chembl_id,
            source=source,
            identifier=identifier,
            selected=primary,
            discarded=candidates[1:],
        )
    return primary


def resolve_pubchem_cid(
    row: Mapping[str, object],
    cache: MutableMapping[str, str | None],
    cfg: PubChemCfg,
) -> str | None:
    """Resolve the preferred PubChem CID for ``row`` using multiple identifiers."""

    chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
    if chembl_id is not None:
        cached = cache.get(chembl_id, _CACHE_MISSING)  # type: ignore[arg-type]
        if cached is not _CACHE_MISSING:
            return cache[chembl_id]

    def _store(resolved: str | None) -> str | None:
        if chembl_id is not None:
            cache[chembl_id] = resolved
        return resolved

    inchikey = _normalise_identifier(
        row.get("standard_inchi_key"), uppercase=True
    )
    if inchikey:
        cid = _select_primary_cid(
            pl.get_cid_from_inchikey(inchikey, cfg),
            chembl_id=chembl_id,
            identifier=inchikey,
            source="standard_inchi_key",
        )
        if cid:
            return _store(cid)

    inchi = _normalise_identifier(row.get("standard_inchi"))
    if inchi:
        cid = _select_primary_cid(
            pl.get_cid_from_inchi(inchi, cfg),
            chembl_id=chembl_id,
            identifier=inchi,
            source="standard_inchi",
        )
        if cid:
            return _store(cid)

    pref_name = _normalise_identifier(row.get("pref_name"))
    if pref_name:
        cid = _select_primary_cid(
            pl.get_cid(pref_name, cfg),
            chembl_id=chembl_id,
            identifier=pref_name,
            source="pref_name_exact",
        )
        if cid:
            return _store(cid)
        cid = _select_primary_cid(
            pl.get_all_cid(pref_name, cfg),
            chembl_id=chembl_id,
            identifier=pref_name,
            source="pref_name_partial",
        )
        if cid:
            return _store(cid)

    smiles = _normalise_identifier(row.get("canonical_smiles"))
    if smiles:
        cid = _select_primary_cid(
            pl.get_cid_from_smiles(smiles, cfg),
            chembl_id=chembl_id,
            identifier=smiles,
            source="canonical_smiles",
        )
        if cid:
            return _store(cid)

    return _store(None)
def attach_parent_molecule_ids(
    df: pd.DataFrame,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None,
    catalog: Mapping[str, str] | None = None,
    source: str | None = None,
) -> tuple[pd.DataFrame, ParentLookupStats]:
    """Attach parent molecule identifiers using the ChEMBL catalogue.

    The function operates on the original ``df`` content and can enrich frames
    that already contain parent identifiers as well as frames without the
    ``catalog_cfg.parent_field`` column.
    """

    if "parant_molecule_id" in df.columns:
        raise ValueError(
            "Input frame contains unexpected column 'parant_molecule_id'."
        )

    result = df.copy()

    if result.empty:
        if catalog_cfg.parent_field not in result.columns:
            result[catalog_cfg.parent_field] = pd.Series(
                pd.NA, index=result.index, dtype="string"
            )
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=0,
            unique=0,
            attached=0,
        )
        return result, stats

    child_column = catalog_cfg.child_field
    parent_column = catalog_cfg.parent_field

    if child_column not in result.columns:
        logger.warning("parent_lookup_missing_child_column", column=child_column)
        result[parent_column] = pd.Series(pd.NA, index=result.index, dtype="string")
        stats = ParentLookupStats(
            source=PARENT_LOOKUP_SOURCE_SKIPPED,
            missing=len(result),
            unique=0,
            attached=0,
        )
        return result, stats

    source_resolved = source
    normalised_child = _normalise_chembl_ids(result[child_column])
    unique_children = normalised_child[normalised_child != ""].unique()
    catalog_data: MutableMapping[str, str] | None

    if catalog is not None:
        base_view = {
            key: catalog[key]
            for key in unique_children
            if key in catalog
        }
        catalog_data = ChainMap({}, base_view)
    else:
        catalog_data = None
    used_partial_cache = False

    if catalog_data is None:
        cache_before = _cache_state(catalog_cfg.cache_path)
        queried = query_parent_catalog(unique_children, catalog_cfg)
        cache_after = _cache_state(catalog_cfg.cache_path)
        if queried or catalog_cfg.sqlite_path.is_file():
            catalog_data = dict(queried)
            source_resolved = _resolve_parent_source(cache_before, cache_after)
            used_partial_cache = True
        else:
            catalog_data = load_parent_catalog(
                client=client,
                api_cfg=api_cfg,
                catalog_cfg=catalog_cfg,
                timeout=timeout,
            )
            cache_after = _cache_state(catalog_cfg.cache_path)
            source_resolved = _resolve_parent_source(cache_before, cache_after)
    else:
        if source_resolved is None:
            cache_exists = catalog_cfg.cache_path.is_file()
            source_resolved = (
                PARENT_LOOKUP_SOURCE_CACHE
                if cache_exists
                else PARENT_LOOKUP_SOURCE_REMOTE
            )

    parent_map = {
        key: catalog_data[key]
        for key in unique_children
        if key in catalog_data
    }
    missing_ids = [key for key in unique_children if key not in parent_map]
    fetched_remote = False

    if missing_ids and catalog is None:
        try:
            fetched = molecule_catalog.fetch_parent_catalog_for(
                missing_ids,
                client=client,
                api_cfg=api_cfg,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("parent_lookup_partial_fetch_failed", error=str(exc))
            fetched = {}
        else:
            if fetched:
                fetched_remote = True
        if fetched:
            catalog_data.update(fetched)
            parent_map.update(fetched)
            if catalog is None:
                if used_partial_cache:
                    update_parent_catalog_cache(fetched, catalog_cfg)
                else:
                    write_parent_catalog_cache(catalog_data, catalog_cfg)

    parent_series = normalised_child.map(parent_map).astype("string")

    if parent_column in result.columns:
        existing_parent = result[parent_column].astype("string")
    else:
        existing_parent = pd.Series(pd.NA, index=result.index, dtype="string")

    combined_parent = existing_parent.copy()
    update_mask = combined_parent.isna() | combined_parent.eq("")
    combined_parent.loc[update_mask] = parent_series.loc[update_mask]
    result[parent_column] = combined_parent.astype("string")

    missing = int(combined_parent.isna().sum())
    attached = len(result) - missing

    if source_resolved is None:
        source_resolved = PARENT_LOOKUP_SOURCE_REMOTE

    stats = ParentLookupStats(
        source=(
            PARENT_LOOKUP_SOURCE_REMOTE
            if fetched_remote
            else source_resolved
        ),
        missing=missing,
        unique=int(len(unique_children)),
        attached=int(attached),
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
    )

    return result, stats


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    cache_path: Path | None = None,
) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    For each canonical SMILES string in ``df``, the function looks up the
    corresponding PubChem CID and basic chemical properties. The PubChem
    fields are appended to the input frame. If a SMILES string cannot be
    resolved, empty values are inserted.

    Parameters
    ----------
    df:
        Data frame returned by :func:`library.chembl_library.get_testitem`.
    cfg:
        PubChem configuration options.

    Returns
    -------
    pandas.DataFrame
        ``df`` with additional PubChem columns.

    """
    if df.empty:
        return df

    result = df.reset_index(drop=True).copy()

    for column in PUBCHEM_COLUMNS:
        if column not in result.columns:
            result[column] = pd.Series(pd.NA, index=result.index, dtype="string")
        else:
            result[column] = result[column].astype("string")

    prefer_local = getattr(cfg, "prefer_local_smiles", False)
    skip_rows: set[int] = set()
    if prefer_local:
        existing_cols = [col for col in PUBCHEM_COLUMNS if col in result.columns]
        if existing_cols:
            normalised = pd.DataFrame(
                {col: result[col].astype("string").replace("", pd.NA) for col in existing_cols}
            )
            complete_mask = normalised.notna().all(axis=1)
            skip_rows = set(normalised[complete_mask].index)

    cache: dict[str, str | None]
    if cache_path is not None:
        cache = _load_pubchem_cid_cache(cache_path)
    else:
        cache = {}

    cache_dirty = False

    if "molecule_chembl_id" in result.columns:
        chembl_ids = result["molecule_chembl_id"].astype("string")
        existing_cids = (
            result["pubchem_cid"].astype("string")
            if "pubchem_cid" in result.columns
            else pd.Series(dtype="string")
        )
        for idx, chembl_value in chembl_ids.items():
            chembl_id = _normalise_identifier(chembl_value, uppercase=True)
            if not chembl_id:
                continue
            existing = (
                _normalise_identifier(existing_cids.iloc[idx])
                if idx < len(existing_cids)
                else None
            )
            if existing is not None and cache.get(chembl_id) != existing:
                cache[chembl_id] = existing
                cache_dirty = True

    indices = [idx for idx in result.index if idx not in skip_rows]
    total = len(indices)
    if total:
        logger.info("pubchem_start", total=total)
    else:
        logger.info("pubchem_no_smiles")

    def _value_or_na(value: str | None) -> object:
        return value if value not in (None, "") else pd.NA

    properties_cache: dict[str, dict[str, object]] = {}
    existing_columns = set(result.columns)

    for position, idx in enumerate(indices, start=1):
        logger.info("pubchem_progress", current=position, total=total)
        row = result.loc[idx]
        chembl_id = _normalise_identifier(row.get("molecule_chembl_id"), uppercase=True)
        previous = cache.get(chembl_id, _CACHE_MISSING) if chembl_id else _CACHE_MISSING
        cid = resolve_pubchem_cid(row, cache, cfg)
        if chembl_id:
            after = cache.get(chembl_id, _CACHE_MISSING)
            if after != previous:
                cache_dirty = True
        if not cid:
            continue
        record = properties_cache.get(cid)
        if record is None:
            props = pl.get_properties(cid, cfg)
            record = {
                "pubchem_cid": cid,
                "pubchem_iupac_name": _value_or_na(props.IUPACName),
                "pubchem_molecular_formula": _value_or_na(props.MolecularFormula),
                "pubchem_isomeric_smiles": _value_or_na(props.iSMILES),
                "pubchem_canonical_smiles": _value_or_na(props.cSMILES),
                "pubchem_inchi": _value_or_na(props.InChI),
                "pubchem_inchikey": _value_or_na(props.InChIKey),
            }
            properties_cache[cid] = record
        for column, value in record.items():
            if prefer_local and column in existing_columns:
                existing_value = result.at[idx, column]
                if not (pd.isna(existing_value) or existing_value == ""):
                    continue
            result.at[idx, column] = value

    for column in PUBCHEM_COLUMNS:
        result[column] = result[column].astype("string")

    if cache_path is not None and cache_dirty:
        _dump_pubchem_cid_cache(cache, cache_path)

    return result


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    limit = cfg.testitem.limit
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="testitem.limit",
            limit=limit,
        )
        return 1

    # Initialise HTTP sessions for downstream HTTP calls
    pl.init_session(cfg.api, cfg.retry)
    # Initialise HTTP session for subsequent ChEMBL requests
    parent_stats = ParentLookupStats(
        source=PARENT_LOOKUP_SOURCE_SKIPPED,
        missing=0,
        unique=0,
        attached=0,
    )
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.testitem.column, cfg=cfg.io
            )
            if limit is not None:
                ids = list(islice(ids_iter, limit))
                logger.info("process_limit", limit=len(ids))
            else:
                ids = list(ids_iter)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        logger.info("identifiers_retrieved", count=len(ids))
        logger.info("chembl_fetch_start", batch_size=cfg.testitem.batch_size)

        try:
            df = cl.get_testitem(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "testitem_fetch_failed",
                error=str(exc),
                batch_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
            )
            return 1
        logger.info("chembl_fetch_done", rows=len(df))
        parent_column = cfg.molecule_catalog.parent_field
        child_column = cfg.molecule_catalog.child_field

        if child_column in df.columns:
            normalised_ids = _normalise_chembl_ids(df[child_column])
        else:
            normalised_ids = pd.Series("", index=df.index, dtype="string")

        if parent_column in df.columns:
            existing_parent = _normalise_chembl_ids(df[parent_column])
        else:
            existing_parent = pd.Series("", index=df.index, dtype="string")

        need_lookup_mask = (normalised_ids != "") & (existing_parent == "")
        need_lookup = set(normalised_ids[need_lookup_mask])

        cache_before = _cache_state(cfg.molecule_catalog.cache_path)
        cache_after = cache_before
        parent_catalog: dict[str, str] = {}
        parent_catalog_source = PARENT_LOOKUP_SOURCE_SKIPPED

        if need_lookup and cache_before[0]:
            try:
                parent_catalog = query_parent_catalog(
                    need_lookup,
                    catalog_cfg=cfg.molecule_catalog,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error(
                    "parent_catalog_invalid",
                    error=str(exc),
                    path=str(cfg.molecule_catalog.cache_path),
                )
                return 1
            cache_after = _cache_state(cfg.molecule_catalog.cache_path)
            parent_catalog_source = _resolve_parent_source(cache_before, cache_after)
            if parent_catalog:
                need_lookup -= set(parent_catalog)

        fetched_remote = False
        if need_lookup:
            try:
                fetched = molecule_catalog.fetch_parent_catalog_for(
                    need_lookup,
                    client=client,
                    api_cfg=cfg.api,
                    timeout=cfg.testitem.timeout,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.error("parent_lookup_partial_fetch_failed", error=str(exc))
                return 1
            if fetched:
                parent_catalog.update(fetched)
                update_parent_catalog_cache(fetched, cfg.molecule_catalog)
                fetched_remote = True
                parent_catalog_source = PARENT_LOOKUP_SOURCE_REMOTE
        if fetched_remote:
            parent_catalog_source = PARENT_LOOKUP_SOURCE_REMOTE
        logger.info("pubchem_augment_start")
        cache_file = cfg.io.cache_dir / PUBCHEM_CID_CACHE_FILENAME
        df = add_pubchem_data(df, cfg.pubchem, cache_path=cache_file)
        logger.info("pubchem_augment_done")
 
        try:
            ensure_no_parant_column(df)
        except ValueError as exc:
            logger.error(
                "invalid_column",
                column=_TYPO_PARENT_COLUMN,
                error=str(exc),
            )
            return 1
 
        logger.info("parent_lookup_start")
        try:
            df, parent_stats = attach_parent_molecule_ids(
                df,
                client=client,
                api_cfg=cfg.api,
                catalog_cfg=cfg.molecule_catalog,
                timeout=cfg.testitem.timeout,
                catalog=parent_catalog,
                source=parent_catalog_source,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("parent_lookup_failed", error=str(exc))
            return 1
        logger.info(
            "parent_lookup_done",
            source=parent_stats.source,
            unique=parent_stats.unique,
            attached=parent_stats.attached,
            missing=parent_stats.missing,
        )

        enrichment_cfg = cfg.testitem_molecule_enrichment
        if enrichment_cfg.enable:
            logger.info("testitem_enrichment_start")
            try:
                df = testitem_enrichment.enrich(
                    df,
                    cfg=enrichment_cfg,
                    io_cfg=cfg.io,
                )
            except ValueError as exc:
                logger.error("testitem_enrichment_failed", error=str(exc))
                return 1
            logger.info("testitem_enrichment_done")

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_testitems(df)
        df = add_pipeline_metadata(df)
        # Determine column order: schema columns first, followed by
        # additional fields sorted alphabetically.
        schema_cols = list(TestitemsSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        rows_total = len(df)
        exit_code = 0
    required_cols = {
        name for name, col in TestitemsSchema.columns.items() if col.required
    }
    optional_cols = set(TestitemsSchema.columns) - required_cols
    missing_required = required_cols - set(df.columns)
    missing_optional = optional_cols - set(df.columns)
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
        try:
            validation_result = validate_testitems(df, return_result=True)
        except SchemaErrors as exc:
            failure_path = Path(output).with_name(
                f"{Path(output).stem}_failure_cases.csv"
            )
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation_failed",
                failures=len(exc.failure_cases),
                path=str(failure_path),
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
        else:
            df = validation_result.data
            if not validation_result.failure_cases.empty:
                failure_path = Path(output).with_name(
                    f"{Path(output).stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in validation_result.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation_failed",
                    failures=len(validation_result.failure_cases),
                    path=str(failure_path),
                )
                exit_code = 1
    else:
        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = ["molecule_chembl_id"]
        csv_path = io.write_csv(
            df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
            col_order=col_order,
        )
        logger.info("write_done", rows=rows_kept, path=str(csv_path))
    except OSError as exc:
        logger.error(
            "write_fail",
            error=str(exc),
            path=str(output),
        )
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="TestitemsSchema",
    )
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error(
            "quality_report_failed",
            error=str(exc),
            path=str(output),
        )
        return 1
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "testitem.timeout",
                "column": "testitem.column",
                "batch_size": "testitem.batch_size",
                "limit": "testitem.limit",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
