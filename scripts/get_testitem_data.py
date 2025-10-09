"""Command line interface for retrieving ChEMBL test item data.

The module wraps :func:`library.pipelines.testitem.run_testitem_pipeline` while
exposing helpers that tests can import directly. Entry points return numeric
exit codes rather than terminating the interpreter to simplify orchestration.
The :func:`ensure_no_parant_column` helper guards against legacy CSV exports
that still include the misspelled parent identifier column.
"""

from __future__ import annotations

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

import argparse
from collections import ChainMap
from collections.abc import Hashable, Mapping, MutableMapping, Sequence
from pathlib import Path
from typing import NamedTuple, cast

import pandas as pd
import requests

from library import (
    cli,  # noqa: F401 - re-exported for monkeypatching in tests
    io,
)
from library.cli import ConfigMetadata, LoggerConfig
from library.cli import build_parser as base_parser
from library.cli.logging import setup_cli_logging
from library.cli.metadata import prepare_option
from library.cli_utils import run_cli_command
from library.clients import pubchem as pc  # noqa: F401 - patched in tests
from library.common.log import logger
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    MoleculeCatalogCfg,
    PubChemCfg,
)
from library.integration import molecule_catalog
from library.integration import pubchem_library as pl
from library.integration.chembl_client import ChemblClient
from library.pipelines import testitem as pipeline
from library.pipelines.testitem import (
    _DEFAULT_CATALOG_CFG,
    PARENT_LOOKUP_SOURCE_CACHE,
    PARENT_LOOKUP_SOURCE_LOOKUP,
    PARENT_LOOKUP_SOURCE_PARTIAL,
    PARENT_LOOKUP_SOURCE_SKIPPED,
    PARENT_LOOKUP_SOURCE_SYNC,
    ParentLookupStats,
    TestitemPipelineOptions,
    load_parent_catalog,
    query_parent_catalog,
    run_testitem_pipeline,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)
from library.pipelines.testitem import catalog as pipeline_catalog
from library.pipelines.testitem import pubchem as pipeline_pubchem

# ===== Parameters =====

DEFAULT_INPUT_NAME = "testitem.csv"
DEFAULT_OUTPUT_STEM = "testitems"

configure_logger = cli.configure_logger

LoadMoleculeHierarchyLookup = pipeline.LoadMoleculeHierarchyLookup

_normalise_identifier = pipeline_pubchem._normalise_identifier
_pubchem_identifiers = pipeline_pubchem._pubchem_identifiers
_pubchem_resolution_key = pipeline_pubchem._pubchem_resolution_key
_load_pubchem_cid_cache = pipeline_pubchem._load_pubchem_cid_cache
resolve_pubchem_cid = pipeline_pubchem.resolve_pubchem_cid

configure_logger = cli.configure_logger


def _resolve_catalog_load_source(
    before: tuple[bool, float | None], after: tuple[bool, float | None]
) -> str:
    """Determine the lookup source after invoking ``load_parent_catalog``."""

    before_exists, before_mtime = before
    after_exists, after_mtime = after

    if after_exists and before_exists and after_mtime == before_mtime:
        return PARENT_LOOKUP_SOURCE_CACHE
    if after_exists:
        return PARENT_LOOKUP_SOURCE_SYNC
    if before_exists and not after_exists:
        return PARENT_LOOKUP_SOURCE_SYNC
    return PARENT_LOOKUP_SOURCE_SYNC


def _normalise_chembl_ids(series: pd.Series) -> pd.Series:
    """Return ``series`` normalised to upper-case ChEMBL identifiers."""

    normalised = series.astype("string").fillna("").str.strip().str.upper()
    return normalised


def load_molecule_hierarchy_lookup(
    path: Path | None,
    *,
    io_cfg: IoCfg,
    encoding: str | None = None,
    delimiter: str | None = None,
    catalog_cfg: MoleculeCatalogCfg | None = None,
) -> dict[str, str | None]:
    """Return child → parent mapping loaded from a local hierarchy file.

    The returned mapping mirrors
    :func:`library.pipelines.testitem.load_molecule_hierarchy_lookup` where
    values may be ``None`` when no parent is listed in the hierarchy.
    """

    cfg_source = catalog_cfg or _DEFAULT_CATALOG_CFG
    resolved_path_value = path or cfg_source.hierarchy_lookup_path
    if resolved_path_value is None:
        return {}

    resolved_path = Path(resolved_path_value)
    resolved_encoding = (
        encoding
        or getattr(cfg_source, "hierarchy_lookup_encoding", None)
        or io_cfg.csv_encoding
    )
    resolved_delimiter = (
        delimiter
        or getattr(cfg_source, "hierarchy_lookup_delimiter", None)
        or io_cfg.csv_sep
    )

    try:
        raw_lookup = _load_molecule_hierarchy_mapping(
            str(resolved_path),
            resolved_encoding,
            resolved_delimiter,
        )
    except FileNotFoundError:
        logger.info("molecule_hierarchy_lookup_missing", path=str(resolved_path))
        return {}
    except ValueError as exc:
        raise ValueError(f"invalid hierarchy lookup: {exc}") from exc

    lookup = {
        child_id: _normalise_parent_identifier(parent_id, child_id=child_id)
        for child_id, parent_id in raw_lookup.items()
    }
    if not lookup:
        return {}

    attached_rows = sum(1 for value in lookup.values() if value is not None)

    logger.info(
        "molecule_hierarchy_lookup_loaded",
        path=str(resolved_path),
        rows=attached_rows,
    )
    return lookup


class ParentLookupPreparedData(NamedTuple):
    """Container for precomputed parent lookup data."""

    child_ids: pd.Series
    existing_parent_ids: pd.Series
    need_lookup: set[str]


def attach_parent_molecule_ids(
    df: pd.DataFrame,
    *,
    client: ChemblClient,
    api_cfg: ApiCfg,
    catalog_cfg: MoleculeCatalogCfg,
    timeout: float | None,
    catalog: Mapping[str, str] | None = None,
    source: str | None = None,
    precomputed: ParentLookupPreparedData | None = None,
) -> tuple[pd.DataFrame, ParentLookupStats]:
    """Attach parent molecule identifiers using the ChEMBL catalogue."""

    if "parant_molecule_id" in df.columns:
        raise ValueError("Input frame contains unexpected column 'parant_molecule_id'.")

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
            uncovered=0,
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
            uncovered=len(result),
        )
        return result, stats

    source_resolved = source
    if precomputed is not None:
        normalised_child = (
            precomputed.child_ids.reindex(result.index, fill_value="")
            .astype("string")
            .copy()
        )
        existing_parent = (
            precomputed.existing_parent_ids.reindex(result.index)
            .astype("string")
            .copy()
        )
        lookup_mask = normalised_child.isin(precomputed.need_lookup)
    else:
        normalised_child = _normalise_chembl_ids(result[child_column])
        if parent_column in result.columns:
            existing_parent = result[parent_column].astype("string").copy()
        else:
            existing_parent = pd.Series(pd.NA, index=result.index, dtype="string")
        lookup_mask = (normalised_child != "") & (
            existing_parent.isna() | existing_parent.eq("")
        )

    existing_parent_before = existing_parent.copy()

    raw_unique_children = tuple(normalised_child[lookup_mask].unique().tolist())
    unique_children = tuple(
        value for value in raw_unique_children if not pd.isna(value)
    )
    catalog_data: MutableMapping[str, str]
    used_partial_cache = False
    needs_full_sync = False
    partial_fetch_used = False
    full_sync_used = False
    uncovered_children = 0
    filters_exclude_parentless = getattr(
        molecule_catalog,
        "_filters_exclude_parentless",
        lambda cfg: False,
    )
    parentless_filtered = bool(filters_exclude_parentless(catalog_cfg))
    json_cache_exists = catalog_cfg.cache_path.is_file()
    sqlite_exists = catalog_cfg.sqlite_path.is_file()

    from library import testitem_pipeline as pipeline_module

    load_catalog_fn = getattr(
        pipeline_module, "load_parent_catalog", load_parent_catalog
    )
    query_catalog_fn = getattr(
        pipeline_module, "query_parent_catalog", query_parent_catalog
    )
    update_cache_fn = getattr(
        pipeline_module, "update_parent_catalog_cache", update_parent_catalog_cache
    )
    write_cache_fn = getattr(
        pipeline_module, "write_parent_catalog_cache", write_parent_catalog_cache
    )

    if catalog is not None:
        base_view = {key: catalog[key] for key in unique_children if key in catalog}
        catalog_data = ChainMap({}, base_view)
    else:
        catalog_data = {}
        if unique_children:
            queried = query_catalog_fn(unique_children, catalog_cfg)
            if queried:
                catalog_data.update(queried)
                used_partial_cache = True
                if source_resolved is None:
                    source_resolved = PARENT_LOOKUP_SOURCE_CACHE
            else:
                used_partial_cache = sqlite_exists
                if sqlite_exists:
                    if source_resolved is None:
                        source_resolved = PARENT_LOOKUP_SOURCE_CACHE

    parent_map = {
        key: catalog_data[key] for key in unique_children if key in catalog_data
    }
    missing_ids = [key for key in unique_children if key not in parent_map]
    uncovered_children = len(missing_ids)

    fetched: dict[str, str] = {}
    if missing_ids and catalog is None:
        if parentless_filtered:
            logger.info(
                "parent_lookup_partial_fetch_skipped_parentless",
                missing=len(missing_ids),
                identifiers=missing_ids,
            )
        else:
            try:
                fetched = molecule_catalog.fetch_parent_catalog_for(
                    missing_ids,
                    client=client,
                    api_cfg=api_cfg,
                    timeout=timeout,
                    catalog_cfg=catalog_cfg,
                )
            except (requests.RequestException, ValueError) as exc:
                logger.warning("parent_lookup_partial_fetch_failed", error=str(exc))
                fetched = {}
            if fetched:
                partial_fetch_used = True
                catalog_data.update(fetched)
                parent_map.update(fetched)
                missing_ids = [key for key in unique_children if key not in parent_map]
                uncovered_children = len(missing_ids)
                if used_partial_cache:
                    update_cache_fn(fetched, catalog_cfg)
                else:
                    write_cache_fn(catalog_data, catalog_cfg)
                    used_partial_cache = True

    needs_full_sync = catalog is None and uncovered_children > 0

    # Определяем, нужно ли фильтровать parentless элементы
    if catalog is None:
        parentless_filtered = bool(filters_exclude_parentless(catalog_cfg))

    skip_full_sync = (
        missing_ids
        and catalog is None
        and needs_full_sync
        and parentless_filtered
        and not (json_cache_exists or sqlite_exists)
    )

    if skip_full_sync:
        logger.info(
            "parent_lookup_full_sync_skipped_parentless",
            count=len(missing_ids),
            identifiers=missing_ids,
        )
        needs_full_sync = False
        source_resolved = PARENT_LOOKUP_SOURCE_SKIPPED
    elif needs_full_sync and parentless_filtered:
        needs_full_sync = False
        logger.info(
            "parent_lookup_skip_full_sync",
            missing=len(missing_ids),
            parentless_filtered=True,
        )
        if source_resolved is None and not partial_fetch_used:
            source_resolved = (
                PARENT_LOOKUP_SOURCE_CACHE
                if used_partial_cache
                else PARENT_LOOKUP_SOURCE_SKIPPED
            )

    if missing_ids and catalog is None and needs_full_sync:
        cache_before_load = pipeline_catalog._cache_state(catalog_cfg.cache_path)
        cache_after_load = cache_before_load
        try:
            loaded_catalog = load_catalog_fn(
                client=client,
                api_cfg=api_cfg,
                catalog_cfg=catalog_cfg,
                timeout=timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.warning("parent_catalog_full_sync_failed", error=str(exc))
        else:
            catalog_data = loaded_catalog
            cache_after_load = pipeline_catalog._cache_state(catalog_cfg.cache_path)
            full_sync_used = True
            source_resolved = _resolve_catalog_load_source(
                cache_before_load, cache_after_load
            )
            if partial_fetch_used:
                catalog_data.update(fetched)
            parent_map = {
                key: catalog_data.get(key, parent_map.get(key, ""))
                for key in unique_children
                if key in catalog_data or key in parent_map
            }
            missing_ids = [key for key in unique_children if key not in parent_map]
            uncovered_children = len(missing_ids)

    if missing_ids:
        log_missing = logger.warning
        if skip_full_sync and parentless_filtered:
            log_missing = logger.info
        log_missing(
            "parent_lookup_missing_parents",
            count=len(missing_ids),
            identifiers=missing_ids,
        )

    refreshed_parent = normalised_child.map(parent_map).astype("string")
    refreshed_normalised = refreshed_parent.fillna("").astype("string")

    combined_parent = existing_parent.astype("string").copy()

    update_mask = combined_parent.isna() | combined_parent.eq("")
    combined_parent.loc[update_mask] = refreshed_parent.loc[update_mask]

    if getattr(catalog_cfg, "force_refresh_existing", False):
        existing_normalised = combined_parent.fillna("").astype("string")
        mismatch_mask = existing_normalised != refreshed_normalised
        if mismatch_mask.any():
            combined_parent.loc[mismatch_mask] = refreshed_parent.loc[mismatch_mask]
        existing_normalised = combined_parent.fillna("").astype("string")

    result[parent_column] = combined_parent.astype("string")

    final_normalised = combined_parent.fillna("").astype("string")
    previous_normalised = existing_parent_before.fillna("").astype("string")

    fallback_mask = (final_normalised != "") & (
        (previous_normalised == "") | (previous_normalised != final_normalised)
    )
    fallback_attached = int(fallback_mask.sum())
    no_parent_count = int((final_normalised == "").sum())
    attached = int((final_normalised != "").sum())
    missing = no_parent_count

    final_source = source_resolved
    if catalog is not None and not missing_ids:
        if final_source in (
            None,
            PARENT_LOOKUP_SOURCE_CACHE,
            PARENT_LOOKUP_SOURCE_SKIPPED,
        ):
            final_source = PARENT_LOOKUP_SOURCE_LOOKUP
    if full_sync_used:
        final_source = PARENT_LOOKUP_SOURCE_SYNC
    elif partial_fetch_used:
        final_source = PARENT_LOOKUP_SOURCE_PARTIAL
    elif final_source is None:
        final_source = (
            PARENT_LOOKUP_SOURCE_CACHE
            if used_partial_cache or catalog is not None
            else PARENT_LOOKUP_SOURCE_SYNC
        )

    stats = ParentLookupStats(
        source=final_source,
        missing=int(missing),
        unique=int(len(unique_children)),
        attached=int(attached),
        uncovered=int(uncovered_children),
        failed_ids=tuple(missing_ids),
        fallback_attached=int(fallback_attached),
        no_parent=int(no_parent_count),
    )

    logger.info(
        "parent_lookup_progress",
        source=stats.source,
        unique=stats.unique,
        attached=stats.attached,
        missing=stats.missing,
        uncovered=stats.uncovered,
        hierarchy_attached=stats.hierarchy_attached,
        fallback_attached=stats.fallback_attached,
        no_parent=stats.no_parent,
    )

    return result, stats


_normalise_parent_identifier = pipeline_catalog._normalise_parent_identifier
_load_molecule_hierarchy_mapping = pipeline_catalog._load_molecule_hierarchy_mapping


def add_pubchem_data(
    df: pd.DataFrame,
    cfg: PubChemCfg,
    *,
    client: ChemblClient | None = None,
    api_cfg: ApiCfg | None = None,
    timeout: float | None = None,
    cid_cache: MutableMapping[str, str | None] | None = None,
    resolution_cache: MutableMapping[Hashable, pl.PubChemResolution] | None = None,
    parent_record_cache: MutableMapping[str, pd.Series | None] | None = None,
    testitem_fields: Sequence[str] | None = None,
    request_limit: int = 1000,
) -> pd.DataFrame:
    """Augment ChEMBL records with PubChem information.

    Delegates to :func:`library.pipelines.testitem.add_pubchem_data` while
    relaxing the ``resolution_cache`` type to align with
    :func:`library.integration.pubchem_library.resolve_pubchem_record`.
    """

    return pipeline.add_pubchem_data(
        df,
        cfg,
        client=client,
        api_cfg=api_cfg,
        timeout=timeout,
        cid_cache=cid_cache,
        resolution_cache=cast(
            MutableMapping[tuple[str | None, ...], pl.PubChemResolution] | None,
            resolution_cache,
        ),
        parent_record_cache=parent_record_cache,
        testitem_fields=testitem_fields,
        request_limit=request_limit,
    )


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Invoke the reusable test item pipeline with CLI parameters.

    Parameters
    ----------
    cfg : Config
        Application configuration containing ChEMBL, PubChem and IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that identifier loading,
        network requests or CSV export failed inside the test item pipeline.
    """

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
        else:
            output_path = None
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
    if output_path is not None:
        args.output_csv = output_path
    options = TestitemPipelineOptions(
        input_csv=Path(args.input_csv),
        output_csv=output_path,
        limit=getattr(args, "limit", None),
        offset=getattr(args, "offset", None),
    )
    return run_testitem_pipeline(cfg, options)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the test item pipeline handling ``--skip-existing`` semantics."""

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            output_path = Path(legacy_output)
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
    else:
        output_path = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
    args.output_csv = output_path
    if args.skip_existing and output_path.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(output_path))
        return 0
    metadata_obj = getattr(args, "_config_metadata", None)
    if not isinstance(metadata_obj, ConfigMetadata):
        metadata_obj = None
    output_source = "cli" if getattr(args, "final_out", None) else "derived"
    limit_value = getattr(args, "limit", None)
    if limit_value is None:
        limit_value = getattr(cfg.testitem, "limit", None)
    offset_value = getattr(args, "offset", getattr(cfg.testitem, "offset", None))
    logger.info(
        "testitem_pipeline_start",
        input=prepare_option(
            metadata_obj, value=str(args.input_csv), default_source="cli"
        ),
        output=prepare_option(
            metadata_obj,
            value=str(output_path),
            default_source=output_source,
        ),
        limit=prepare_option(
            metadata_obj,
            argument="limit",
            path="sources.chembl.pipelines.testitem.limit",
            value=limit_value,
        ),
        offset=prepare_option(
            metadata_obj,
            argument="offset",
            path="sources.chembl.pipelines.testitem.offset",
            value=offset_value,
        ),
        batch_size=prepare_option(
            metadata_obj,
            argument="batch_size",
            path="sources.chembl.pipelines.testitem.batch_size",
            value=getattr(cfg.testitem, "batch_size", None),
        ),
        timeout=prepare_option(
            metadata_obj,
            argument="timeout",
            path="sources.chembl.pipelines.testitem.timeout",
            value=getattr(cfg.testitem, "timeout", None),
        ),
        column=prepare_option(
            metadata_obj,
            argument="column",
            path="sources.chembl.pipelines.testitem.column",
            value=getattr(cfg.testitem, "column", None),
        ),
    )
    exit_code = run_chembl(cfg, args)
    if exit_code == 0:
        logger.info("testitem_pipeline_done", output=str(output_path))
    else:
        logger.error(
            "testitem_pipeline_failed",
            output=str(output_path),
            exit_code=exit_code,
        )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser configured with the common CLI options and the associated logging
        configuration used by :func:`main`.
    """

    parser, log_cfg = base_parser(
        "ChEMBL and PubChem compound data utilities",
        column="molecule_chembl_id",
        chunk_size=1000,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
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
        help=("Maximum number of identifiers to process; use 0 to skip processing"),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the values from
        :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the pipeline finishes successfully, non-zero otherwise.

    Raises
    ------
    SystemExit
        Raised when invalid command-line options are provided.
    """

    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)

    cli.prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )

    if args.limit == 0:
        logger.info("pipeline_skip_limit", limit=args.limit)
        return 0

    if args.limit is not None and args.limit < 0:
        parser.error("--limit must be zero or a positive integer")
    if args.offset < 0:
        parser.error("--offset must be zero or a positive integer")

    mapping = {
        "timeout": "testitem.timeout",
        "column": "testitem.column",
        "batch_size": "testitem.batch_size",
        "limit": "testitem.limit",
        "offset": "testitem.offset",
    }
    with setup_cli_logging(
        Path(__file__).with_suffix("").name, log_cfg, getattr(args, "date", None)
    ) as logging_ctx:
        exit_code = run_cli_command(
            args=args,
            parser=parser,
            log_cfg=logging_ctx.log_cfg,
            mapping=mapping,
            run=run,
            logger=logger,
        )
    configure_logger(log_cfg)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
