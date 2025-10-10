"""Enrich test item data with salt information and molecule catalog flags."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ... import io
from ...common.log import logger
from ...config import IoCfg, TestitemMoleculeEnrichmentCfg

_FLAG_COLUMNS: tuple[str, ...] = ("natural_product", "prodrug", "polymer_flag")
_TRUE_VALUES: frozenset[str] = frozenset({"1", "true", "t", "yes", "y"})
_FALSE_VALUES: frozenset[str] = frozenset({"0", "false", "f", "no", "n"})


@dataclass(frozen=True)
class _SourceFrames:
    hierarchy: pd.DataFrame
    catalog: pd.DataFrame


def _read_source(path: Path, *, io_cfg: IoCfg) -> pd.DataFrame:
    """Read a CSV file using the shared I/O configuration."""

    try:
        return io.read_csv(path, cfg=io_cfg, dtype=pd.StringDtype())
    except io.CsvReadError as exc:
        logger.error(
            "testitem_enrichment_source_read_failed",
            path=str(path),
            error=str(exc.original_error),
        )
        raise


def _normalise_ids(series: pd.Series) -> pd.Series:
    """Return ``series`` with whitespace trimmed and converted to upper case."""

    normalised = series.fillna("").astype("string").str.strip().str.upper()
    return normalised.mask(normalised == "", other=pd.NA)


def _coerce_flag_values(series: pd.Series, column: str) -> tuple[pd.Series, set[str]]:
    """Coerce string flags to pandas nullable booleans."""

    unknown: set[str] = set()
    result: list[object] = []
    for value in series.astype(object).tolist():
        if value is None or pd.isna(value):
            result.append(pd.NA)
            continue
        text = str(value).strip()
        if not text:
            result.append(pd.NA)
            continue
        text_lower = text.lower()
        if text_lower in _TRUE_VALUES:
            result.append(True)
        elif text_lower in _FALSE_VALUES:
            result.append(False)
        else:
            result.append(pd.NA)
            unknown.add(text_lower)
    return pd.Series(result, index=series.index, dtype="boolean"), unknown


def _summarise_identifiers(
    values: Sequence[str], *, limit: int = 20
) -> tuple[list[str], bool]:
    """Return a truncated list of identifiers and whether truncation occurred."""

    identifiers = list(values)
    truncated = len(identifiers) > limit
    if truncated:
        identifiers = identifiers[:limit]
    return identifiers, truncated


def _load_sources(
    cfg: TestitemMoleculeEnrichmentCfg, *, io_cfg: IoCfg
) -> _SourceFrames:
    """Load hierarchy and catalog frames according to configuration."""

    hierarchy_required = {"molecule_chembl_id", "parent_molecule_chembl_id"}
    catalog_required = {"molecule_chembl_id", *_FLAG_COLUMNS}

    hierarchy = pd.DataFrame(columns=list(hierarchy_required))
    catalog = pd.DataFrame(columns=list(catalog_required))

    try:
        hierarchy = _read_source(cfg.sources.molecule_hierarchy_path, io_cfg=io_cfg)
    except FileNotFoundError:
        logger.warning(
            "testitem_enrichment_missing_hierarchy",
            path=str(cfg.sources.molecule_hierarchy_path),
        )
    else:
        missing = hierarchy_required - set(hierarchy.columns)
        if missing:
            raise ValueError(
                "missing columns in hierarchy file: " + ", ".join(sorted(missing))
            )
        hierarchy = hierarchy.loc[:, list(hierarchy_required)].copy()

    try:
        catalog = _read_source(cfg.sources.molecule_catalog_path, io_cfg=io_cfg)
    except FileNotFoundError:
        logger.warning(
            "testitem_enrichment_missing_catalog",
            path=str(cfg.sources.molecule_catalog_path),
        )
    else:
        missing = catalog_required - set(catalog.columns)
        if missing:
            raise ValueError(
                "missing columns in molecule catalog: " + ", ".join(sorted(missing))
            )
        catalog = catalog.loc[:, list(catalog_required)].copy()

    return _SourceFrames(hierarchy=hierarchy, catalog=catalog)


def _prepare_catalog(
    catalog: pd.DataFrame,
    cfg: TestitemMoleculeEnrichmentCfg,
) -> tuple[pd.DataFrame, dict[str, set[str]]]:
    """Return normalised catalog frame and unknown flag values."""

    unknown_by_column: dict[str, set[str]] = {column: set() for column in _FLAG_COLUMNS}

    if catalog.empty:
        return catalog, unknown_by_column

    catalog["molecule_chembl_id"] = _normalise_ids(catalog["molecule_chembl_id"])
    catalog = catalog.dropna(subset=["molecule_chembl_id"]).drop_duplicates(
        subset=["molecule_chembl_id"], keep="first"
    )

    if cfg.flags.coerce_to_bool:
        for column in _FLAG_COLUMNS:
            catalog[column], unknown_values = _coerce_flag_values(
                catalog[column], column
            )
            unknown_by_column[column] = unknown_values
    else:
        for column in _FLAG_COLUMNS:
            catalog[column] = catalog[column].astype("string")

    return catalog.set_index("molecule_chembl_id"), unknown_by_column


def _map_flags(ids: pd.Series, mapping: pd.Series) -> pd.Series:
    """Map identifiers to boolean flag values using the provided mapping."""

    values: list[object] = []
    lookup = mapping.to_dict()
    for value in ids.astype(object).tolist():
        if value is None:
            values.append(pd.NA)
        else:
            values.append(lookup.get(str(value), pd.NA))
    return pd.Series(values, index=ids.index, dtype="boolean")


def enrich(
    df: pd.DataFrame,
    *,
    cfg: TestitemMoleculeEnrichmentCfg,
    io_cfg: IoCfg,
) -> pd.DataFrame:
    """Attach salt identifiers and molecule flags to ``df``."""

    if df.empty:
        df = df.copy()
        df["salt_chembl_id"] = pd.Series(dtype="string")
        for column in _FLAG_COLUMNS:
            df[column] = pd.Series(dtype="boolean")
        return df

    frames = _load_sources(cfg, io_cfg=io_cfg)

    result = df.copy()
    if "molecule_chembl_id" in result.columns:
        child_ids = _normalise_ids(result["molecule_chembl_id"])
    else:
        child_ids = pd.Series(pd.NA, index=result.index, dtype="string")

    if "parent_molecule_chembl_id" in result.columns:
        parent_ids = _normalise_ids(result["parent_molecule_chembl_id"])
    else:
        parent_ids = pd.Series(pd.NA, index=result.index, dtype="string")
        result["parent_molecule_chembl_id"] = parent_ids

    if not frames.hierarchy.empty:
        hierarchy = frames.hierarchy.copy()
        hierarchy["molecule_chembl_id"] = _normalise_ids(
            hierarchy["molecule_chembl_id"]
        )
        hierarchy["parent_molecule_chembl_id"] = _normalise_ids(
            hierarchy["parent_molecule_chembl_id"]
        )
        hierarchy = hierarchy.dropna(subset=["molecule_chembl_id"])
        hierarchy = hierarchy.drop_duplicates(
            subset=["molecule_chembl_id"], keep="first"
        ).set_index("molecule_chembl_id")
        parent_from_hierarchy = child_ids.map(
            hierarchy["parent_molecule_chembl_id"].to_dict()
        )
        parent_ids = parent_ids.fillna(parent_from_hierarchy)

    parent_ids = parent_ids.astype("string")
    result["parent_molecule_chembl_id"] = parent_ids

    salt_mask = parent_ids.notna() & (parent_ids != child_ids)
    salt_series = pd.Series(pd.NA, index=result.index, dtype="string")
    salt_series[salt_mask] = child_ids[salt_mask]
    if not cfg.output.salt_as_null_when_absent:
        salt_series = salt_series.fillna("-")
    result["salt_chembl_id"] = salt_series

    catalog, unknown_flags = _prepare_catalog(frames.catalog, cfg)

    if cfg.logging.warn_missing_molecule and not catalog.empty:
        known_ids = set(catalog.index)
        missing_children = sorted(
            {value for value in child_ids.dropna() if value not in known_ids}
        )
        missing_parent_values = {
            value for value in parent_ids.dropna() if value not in known_ids
        }
        missing_parentless_children: set[str] = set()
        for child_value, parent_value in zip(child_ids, parent_ids, strict=False):
            if pd.isna(parent_value) and pd.notna(child_value):
                child_text = str(child_value)
                if child_text not in known_ids:
                    missing_parentless_children.add(child_text)
        missing_parents = sorted(missing_parent_values | missing_parentless_children)
        if missing_children:
            child_identifiers, child_truncated = _summarise_identifiers(
                missing_children
            )
            logger.warning(
                "testitem_enrichment_missing_child_flags",
                count=len(missing_children),
                identifiers=child_identifiers,
                truncated=child_truncated,
            )
        if missing_parents:
            parent_identifiers, parent_truncated = _summarise_identifiers(
                missing_parents
            )
            logger.warning(
                "testitem_enrichment_missing_parent_flags",
                count=len(missing_parents),
                identifiers=parent_identifiers,
                truncated=parent_truncated,
            )

    for column, values in unknown_flags.items():
        if values and cfg.flags.coerce_to_bool:
            logger.warning(
                "testitem_enrichment_unknown_flag_values",
                column=column,
                values=sorted(values),
            )

    if catalog.empty:
        for column in _FLAG_COLUMNS:
            result[column] = pd.Series(pd.NA, index=result.index, dtype="boolean")
        return result

    for column in _FLAG_COLUMNS:
        mapping = catalog[column]
        if cfg.flags.coerce_to_bool:
            child_flag = _map_flags(child_ids, mapping)
            parent_flag = _map_flags(parent_ids, mapping)
        else:
            lookup = mapping.to_dict()
            child_flag = child_ids.map(lookup).astype("string")
            parent_flag = parent_ids.map(lookup).astype("string")

        if cfg.logging.warn_inconsistent_flags:
            inconsistent = (
                child_flag.notna() & parent_flag.notna() & (child_flag != parent_flag)
            )
            if inconsistent.any():
                logger.warning(
                    "testitem_enrichment_inconsistent_flag",
                    column=column,
                    count=int(inconsistent.sum()),
                )

        if cfg.flags.parent_fallback:
            missing_mask = child_flag.isna()
            child_flag[missing_mask] = parent_flag[missing_mask]

        if cfg.flags.coerce_to_bool:
            result[column] = child_flag.astype("boolean")
        else:
            result[column] = child_flag

    return result


__all__ = ["enrich"]
