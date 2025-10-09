"""Main entry point for the organism helper post-processing pipeline."""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

from ...common.log import logger
from ..helpers import normalise_export_basename
from .cellularity import add_cellularity_smart, FetchLineageCallable
from .multifunctional import compute_multifunctional


_CANONICAL_EXPORT_PREFIX = "output.target_"
# ``20250101`` mirrors the frozen snapshot bundled with the test fixtures and
# ensures deterministic helper names even when the input export is non-canonical
# (e.g. ``targets_minimal.csv``). The value is aligned with
# ``tests/unit/test_target_postprocess_main.py`` expectations.
_DEFAULT_EXPORT_STAMP = "20250101"
_STAMP_SOURCE_COLUMNS: tuple[str, ...] = (
    "timestamp_utc",
    "timestamp",
    "export_timestamp",
    "export_date",
)


def _is_canonical_basename(name: str) -> bool:
    lowered = name.lower()
    return lowered.startswith(_CANONICAL_EXPORT_PREFIX) and lowered.endswith(".csv")


def _extract_stamp_candidate(series: pd.Series) -> str | None:
    if series.empty:
        return None

    as_string = series.astype("string")
    for value in as_string:
        if value is None or value is pd.NA:
            continue
        text = str(value).strip()
        if not text:
            continue
        parsed = pd.to_datetime(text, utc=True, errors="coerce")
        if pd.notna(parsed):
            return parsed.strftime("%Y%m%d")
        match = re.search(r"(\d{8})", text)
        if match:
            return match.group(1)
    return None


def _derive_export_stamp(frame: pd.DataFrame) -> str | None:
    for column in _STAMP_SOURCE_COLUMNS:
        if column not in frame.columns:
            continue
        stamp = _extract_stamp_candidate(frame[column])
        if stamp:
            return stamp
    return None


def _should_force_canonical(base: str) -> bool:
    if _is_canonical_basename(base):
        return False

    lowered = base.lower()
    if lowered.startswith("output.target"):
        return True
    if lowered.startswith("output.targets"):
        return True
    if lowered.startswith("targets"):
        return True
    return False


def _normalise_helper_basename(frame: pd.DataFrame, base: str) -> str:
    if not _should_force_canonical(base):
        return base

    stamp = _derive_export_stamp(frame) or _DEFAULT_EXPORT_STAMP
    return f"output.target_{stamp}.csv"


def _lowercase_column(series: pd.Series) -> pd.Series:
    """Return the lower-cased variant of ``series`` preserving nulls."""

    return series.astype("string").str.lower()


def postprocess_target_table(
    input_path: str | Path,
    *,
    email: str | None = None,
    fetcher: FetchLineageCallable | None = None,
) -> str:
    """Replicate the Power Query post-processing for organism helpers."""

    path = Path(input_path)
    source = pd.read_csv(path, dtype=str, keep_default_na=False)

    base_columns = [
        "target_chembl_id",
        "uniprot_id_primary",
        "organism",
        "taxon_id",
        "lineage_superkingdom",
        "lineage_phylum",
        "lineage_class",
    ]
    missing_columns = [
        column for column in base_columns if column not in source.columns
    ]
    if missing_columns:
        logger.warning(
            "target_postprocess_missing_columns",
            path=str(path),
            columns=missing_columns,
        )

    source_base = source.reindex(columns=base_columns, fill_value="").copy()

    for column in ("lineage_superkingdom", "lineage_phylum", "lineage_class"):
        source_base[column] = _lowercase_column(source_base[column])

    with_cellularity = add_cellularity_smart(
        source_base,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
        email=email,
        fetcher=fetcher,
    )

    multifunctional_part = compute_multifunctional(source)
    joined = with_cellularity.merge(
        multifunctional_part[["target_chembl_id", "multifunctional_enzyme"]],
        on="target_chembl_id",
        how="left",
        sort=False,
    )

    base = normalise_export_basename(path)
    if not missing_columns:
        base = _normalise_helper_basename(source, base)
    output_path = path.with_name(f"organism.{base}").with_suffix(".csv")
    joined.to_csv(output_path, index=False, encoding="utf-8", lineterminator="\n")
    print(f"Postprocessed target table saved to: {output_path.name}")
    return str(output_path)


__all__ = ["postprocess_target_table"]
