"""Activity post-processing mirroring the legacy Power Query workbook.

The historical ChEMBL acquisition pipeline relied on an Excel workbook that
merged activity measurements with dictionary tables describing citation
significance and target metadata.  The workbook emitted
``extended.output.activity_<stamp>.csv`` files combining boolean flags,
taxonomic annotations, and lookup-driven enrichments.  This module is the
Python port of that workflow.  It reproduces every transformation step so that
the resulting CSVs remain byte-identical with the legacy exports.
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from library.common.log import logger
from library.pipelines.target import organism_classification

from . import helpers

_DEFAULT_SEARCH_DIR = Path("data/output")
_DEFAULT_DICTIONARY_DIR = Path("config/dictionary")
_FILENAME_RE = re.compile(r"output\.activit(?:y|ies)_(\d{8})\.csv\Z")

_REQUIRED_COLUMNS: frozenset[str] = frozenset(
    {
        "activity_chembl_id",
        "salt_chembl_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "bao_endpoint",
        "standard_type",
        "standard_value",
        "log_value",
        "bao_format",
        "compound_key",
        "compound_name",
        "multmol_assay",
        "approx_cited_activity",
        "shuffled_cit",
        "exact_cited_activity",
        "higly_correlated_cit",
        "review_doc",
        "rounded_data_citation",
        "original_activity_approx",
        "original_activity_exact",
        "nstereo",
    }
)

_REQUIRED_COLUMN_DTYPES: Mapping[str, str] = {
    "activity_chembl_id": "string",
    "salt_chembl_id": "string",
    "molecule_chembl_id": "string",
    "target_chembl_id": "string",
    "assay_chembl_id": "string",
    "document_chembl_id": "string",
    "bao_endpoint": "string",
    "standard_type": "string",
    "standard_value": "Float64",
    "log_value": "Float64",
    "bao_format": "string",
    "compound_key": "string",
    "compound_name": "string",
    "multmol_assay": "boolean",
    "approx_cited_activity": "boolean",
    "shuffled_cit": "boolean",
    "exact_cited_activity": "boolean",
    "higly_correlated_cit": "boolean",
    "review_doc": "boolean",
    "rounded_data_citation": "boolean",
    "original_activity_approx": "string",
    "original_activity_exact": "string",
    "nstereo": "Int64",
}

_REQUIRED_COLUMN_FALLBACKS: Mapping[str, Callable[[pd.DataFrame], pd.Series | None]] = {
    "activity_chembl_id": lambda frame: frame.get("activity_id"),
    "salt_chembl_id": lambda frame: frame.get("molecule_chembl_id"),
    "compound_name": lambda frame: frame.get("molecule_pref_name"),
    "log_value": lambda frame: frame.get("pchembl_value"),
}

_GROUP_KEY_COLUMNS: tuple[str, ...] = (
    "salt_chembl_id",
    "molecule_chembl_id",
    "target_chembl_id",
    "assay_chembl_id",
    "standard_type",
)

_FINAL_COLUMN_ORDER: tuple[str, ...] = (
    "activity_chembl_id",
    "saltform_id",
    "molecule_chembl_id",
    "target_chembl_id",
    "assay_chembl_id",
    "document_chembl_id",
    "bao_endpoint",
    "standard_type",
    "standard_value",
    "pA_value",
    "bao_format",
    "compound_key",
    "compound_name",
    "unknown_chirality",
    "multmol_assay",
    "exact_data_citation",
    "higly_correlated_assay",
    "shuffled_assay",
    "review",
    "rounded_data_citation",
    "high_citation_rate",
    "original_activity_approx",
    "original_activity_exact",
    "is_citation",
    "IUPHAR_class",
    "IUPHAR_subclass",
    "unicellular_organism",
    "multifunctional_enzyme",
    "gene_index",
    "taxon_index",
    "target_sort_order",
)

_TARGET_COLUMNS: Sequence[str] = (
    "target_chembl_id",
    "target_sort_order",
    "multifunctional_enzyme",
    "IUPHAR_class",
    "IUPHAR_subclass",
    "genus",
    "superkingdom",
    "phylum",
    "taxon_id",
    "gene_index",
    "taxon_index",
)


class ActivityExtendedError(RuntimeError):
    """Raised when the activity extended post-processing cannot proceed."""


def _current_default_search_dir() -> Path:
    package = sys.modules.get(__name__)
    if package is not None and hasattr(package, "_DEFAULT_SEARCH_DIR"):
        override = getattr(package, "_DEFAULT_SEARCH_DIR")
        if override is not None:
            return Path(override)
    return _DEFAULT_SEARCH_DIR


def _normalised_activity_basename(path: Path) -> tuple[str, str] | None:
    """Return ``(date_token, name)`` for ``path`` when it matches the activity pattern."""

    try:
        candidate = helpers.normalise_export_basename(path)
    except ValueError:
        return None
    match = _FILENAME_RE.match(candidate)
    if match is None:
        return None
    return match.group(1), candidate


def _latest_activity_export(search_dir: Path) -> Path:
    candidates: list[tuple[str, str, Path]] = []
    for path in search_dir.iterdir():
        if not path.is_file():
            continue
        normalised = _normalised_activity_basename(path)
        if normalised is None:
            continue
        candidates.append((normalised[0], normalised[1], path))
    candidates.sort(key=lambda item: (item[0], item[1]))
    if not candidates:
        raise ActivityExtendedError(
            "No activity exports matching 'output.activity_YYYYMMDD.csv' or "
            f"'output.activities_YYYYMMDD.csv' found in {search_dir!s}"
        )
    return candidates[-1][2]


def _resolve_dictionary_root(dictionary_dir: Path | None) -> Path:
    if dictionary_dir is None:
        return _DEFAULT_DICTIONARY_DIR
    return Path(dictionary_dir)


def _load_citation_fraction(dictionary_root: Path) -> pd.DataFrame:
    candidate = dictionary_root / "_activity" / "citation_fraction.csv"
    if not candidate.exists():
        raise ActivityExtendedError(
            "citation_fraction.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
        )
    return pd.read_csv(candidate, dtype={"N": "Int64", "K_min_significant": "Int64"})


def _resolve_targets_path(dictionary_root: Path, override: Path | None) -> Path:
    if override is not None:
        path = Path(override)
        if not path.exists():
            raise ActivityExtendedError(f"targets_type.csv override not found: {path!s}")
        return path

    candidates = [
        dictionary_root / "targets_type.csv",
        dictionary_root / "_target" / "targets_type.csv",
        dictionary_root / "_Target" / "targets_type.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    formatted = " or ".join(f"'{path}'" for path in candidates[:-1]) + f" or '{candidates[-1]}'"
    raise ActivityExtendedError(
        "targets_type.csv not found in the provided dictionary directory. Expected at "
        + formatted
    )


def _load_target_metadata(path: Path) -> pd.DataFrame:
    dtype: Mapping[str, str] = {
        "target_chembl_id": "string",
        "IUPHAR_class": "string",
        "IUPHAR_subclass": "string",
        "gene_index": "string",
        "taxon_index": "string",
        "target_sort_order": "string",
        "multifunctional_enzyme": "string",
        "genus": "string",
        "superkingdom": "string",
        "phylum": "string",
        "taxon_id": "string",
    }
    return pd.read_csv(path, usecols=_TARGET_COLUMNS, dtype=dtype)


def _safe_to_bool(series: pd.Series, column: str) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise ActivityExtendedError(f"column '{column}' has duplicate entries; expected a Series")

    def mapper(value: object) -> object:
        if pd.isna(value):
            return pd.NA
        if isinstance(value, str):
            lowered = value.strip().lower()
        else:
            lowered = value
        if lowered in (True, 1, "1", "true", "t"):
            return True
        if lowered in (False, 0, "0", "false", "f"):
            return False
        raise ValueError(f"invalid boolean value: {value}")

    try:
        mapped = series.map(mapper)
        return mapped.astype("boolean")
    except Exception as exc:  # pragma: no cover - defensive downgrade
        logger.warning("dtype_bool_conversion_failed", column=column, error=str(exc))
        return series.astype("string")


def _safe_to_int(series: pd.Series, column: str) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise ActivityExtendedError(f"column '{column}' has duplicate entries; expected a Series")
    try:
        return pd.to_numeric(series, errors="raise").astype("Int64")
    except Exception as exc:  # pragma: no cover - defensive downgrade
        logger.warning("dtype_int_conversion_failed", column=column, error=str(exc))
        return series.astype("string")


def _prepare_unknown_chirality(frame: pd.DataFrame) -> pd.DataFrame:
    df = frame.copy()
    if "nstereo" in df.columns:
        df["unknown_chirality"] = _safe_to_int(df["nstereo"], "nstereo").ne(1).fillna(True)
        df.drop(columns=["nstereo"], inplace=True)
    else:
        df["unknown_chirality"] = pd.Series(True, index=df.index, dtype="boolean")
    return df


def _apply_multimol_logic(df: pd.DataFrame) -> pd.DataFrame:
    missing = set(_GROUP_KEY_COLUMNS) - set(df.columns)
    if missing:
        raise ActivityExtendedError(
            "activity table missing columns for multimol grouping: " + ", ".join(sorted(missing))
        )
    counts = (
        df.groupby(list(_GROUP_KEY_COLUMNS), dropna=False)
        .size()
        .rename("Count")
        .reset_index()
    )
    merged = df.merge(counts, on=list(_GROUP_KEY_COLUMNS), how="left")
    mask = (
        merged["unknown_chirality"].fillna(True).eq(False)
        & merged["multmol_assay"].isna()
        & (merged["Count"] > 1)
    )
    duplicated_assays = set(merged.loc[mask, "assay_chembl_id"].dropna().astype(str))
    merged["multimol_assay_same"] = merged["assay_chembl_id"].isin(duplicated_assays)

    multmol_series = _safe_to_bool(merged["multmol_assay"], "multmol_assay").fillna(False)
    merged["multmol_assay"] = _safe_to_bool(
        multmol_series | merged["multimol_assay_same"], "multmol_assay"
    )
    merged.drop(columns=["multimol_assay_same", "Count"], inplace=True)
    return merged


def _empty_series(index: pd.Index, dtype: str) -> pd.Series:
    if dtype == "boolean":
        return pd.Series(pd.NA, index=index, dtype="boolean")
    if dtype == "Float64":
        return pd.Series(pd.NA, index=index, dtype="Float64")
    if dtype == "Int64":
        return pd.Series(pd.NA, index=index, dtype="Int64")
    return pd.Series(pd.NA, index=index, dtype=dtype)


def _ensure_required_input_columns(frame: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    df = frame.copy()
    filled: set[str] = set()
    if df.empty:
        for column, dtype in _REQUIRED_COLUMN_DTYPES.items():
            if column not in df.columns:
                df[column] = pd.Series([], dtype=dtype)
                filled.add(column)
        return df, filled

    for column, dtype in _REQUIRED_COLUMN_DTYPES.items():
        if column in df.columns:
            continue
        fallback = _REQUIRED_COLUMN_FALLBACKS.get(column)
        if fallback is not None:
            candidate = fallback(df)
            if candidate is not None:
                aligned = candidate.reindex(df.index)
                try:
                    df[column] = aligned.astype(dtype)
                except (TypeError, ValueError):
                    df[column] = aligned.astype("string")
                filled.add(column)
                continue
        df[column] = _empty_series(df.index, dtype)
        filled.add(column)
    return df, filled


def _rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    renamed = df.rename(
        columns={
            "salt_chembl_id": "saltform_id",
            "approx_cited_activity": "rounded_data_citation",
            "shuffled_cit": "shuffled_assay",
            "exact_cited_activity": "exact_data_citation",
            "higly_correlated_cit": "higly_correlated_assay",
            "review_doc": "review",
            "log_value": "pA_value",
        }
    )
    return renamed.loc[:, ~renamed.columns.duplicated(keep="last")]


def _drop_unused_columns(df: pd.DataFrame) -> pd.DataFrame:
    drop_cols = [
        "cited_document",
        "acts_per_assay_step5",
        "approx_cited_activity_samedoc",
        "error_document",
        "exact_cited_activity_samedoc",
        "higly_correlated_cit_samedoc",
        "standard_inchi_skeleton",
        "standard_inchi_stereo",
        "step1",
        "step2",
        "step3",
        "step4",
        "step5",
        "step6",
        "survives_main_steps",
    ]
    present = [column for column in drop_cols if column in df.columns]
    return df.drop(columns=present, errors="ignore")


def _compute_citation_flags(df: pd.DataFrame) -> pd.DataFrame:
    bool_columns = [
        "exact_data_citation",
        "higly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
    ]
    converted = df.copy()
    for column in bool_columns:
        if column in converted.columns:
            converted[column] = _safe_to_bool(converted[column], column).fillna(False)
        else:
            converted[column] = pd.Series(False, index=converted.index, dtype="boolean")
    converted["is_citation"] = converted[bool_columns].any(axis=1)
    return converted


def _annotate_high_citation(df: pd.DataFrame, dictionary_root: Path) -> pd.DataFrame:
    converted = df.copy()
    counts = (
        converted.groupby("document_chembl_id")["is_citation"]
        .agg(n_citation="sum", n_non_citation=lambda s: (~s).sum())
        .reset_index()
    )
    counts["N"] = counts["n_citation"] + counts["n_non_citation"]
    counts = counts[(counts["n_citation"] > 0) & (counts["n_non_citation"] > 0)]

    citation_fraction = _load_citation_fraction(dictionary_root)
    counts = counts.merge(
        citation_fraction[["N", "K_min_significant"]],
        on="N",
        how="left",
    )
    counts["high_citation_rate"] = counts["K_min_significant"].notna() & (
        counts["n_citation"] >= counts["K_min_significant"]
    )

    converted = converted.merge(
        counts[["document_chembl_id", "high_citation_rate"]],
        on="document_chembl_id",
        how="left",
    )
    converted["high_citation_rate"] = _safe_to_bool(
        converted["high_citation_rate"], "high_citation_rate"
    ).fillna(False)
    return converted


def _merge_target_metadata(
    df: pd.DataFrame,
    *,
    dictionary_root: Path,
    targets_override: Path | None,
) -> pd.DataFrame:
    targets_path = _resolve_targets_path(dictionary_root, targets_override)
    targets = _load_target_metadata(targets_path)
    merged = df.merge(targets, on="target_chembl_id", how="left")
    merged = merged.loc[:, ~merged.columns.duplicated()]

    if "organism_cellularity" not in merged.columns:
        merged["organism_cellularity"] = pd.Series(dtype="string")

    merged["multifunctional_enzyme"] = _safe_to_bool(
        merged["multifunctional_enzyme"], "multifunctional_enzyme"
    )

    enriched = organism_classification.add_cellularity_smart(
        merged,
        genus_col="genus",
        superkingdom_col="superkingdom",
        phylum_col="phylum",
        output_col="organism_cellularity",
    )

    unicellular_labels = {
        organism_classification.TYPE_UNICELLULAR,
        organism_classification.TYPE_VIRAL,
    }
    enriched["unicellular_organism"] = (
        enriched["organism_cellularity"].astype("string").isin(unicellular_labels)
    ).astype("boolean")

    enriched.drop(columns=["genus", "superkingdom", "phylum", "taxon_id", "organism_cellularity"], inplace=True)
    return enriched


def _select_and_cast(df: pd.DataFrame) -> pd.DataFrame:
    missing = set(_FINAL_COLUMN_ORDER) - set(df.columns)
    if missing:
        raise ActivityExtendedError(
            "activity table missing expected columns: " + ", ".join(sorted(missing))
        )
    result = df.loc[:, _FINAL_COLUMN_ORDER].copy()
    bool_columns = [
        "unknown_chirality",
        "multmol_assay",
        "exact_data_citation",
        "higly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
        "high_citation_rate",
        "is_citation",
        "unicellular_organism",
        "multifunctional_enzyme",
    ]
    for column in bool_columns:
        if column in result.columns:
            result[column] = _safe_to_bool(result[column], column)
    return result


def _augment_activity_frame(frame: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    df = frame.copy()
    filled: set[str] = set()

    if "activity_chembl_id" not in df.columns and "activity_id" in df.columns:
        df["activity_chembl_id"] = df["activity_id"].astype("string")
        filled.add("activity_chembl_id")

    if "compound_name" not in df.columns and "molecule_pref_name" in df.columns:
        df["compound_name"] = df["molecule_pref_name"].astype("string")
        filled.add("compound_name")

    compound_key_missing: pd.Series | None = None
    if "compound_key" in df.columns:
        compound_key_missing = df["compound_key"].isna()
    for candidate in ("parent_molecule_chembl_id", "molecule_chembl_id"):
        if candidate not in df.columns:
            continue
        candidate_values = df[candidate].astype("string")
        if "compound_key" not in df.columns:
            df["compound_key"] = candidate_values
            filled.add("compound_key")
            break
        assert compound_key_missing is not None  # for type checkers
        if compound_key_missing.any():
            df.loc[compound_key_missing, "compound_key"] = candidate_values.loc[compound_key_missing]
            filled.add("compound_key")
            compound_key_missing = df["compound_key"].isna()
        if compound_key_missing is not None and not compound_key_missing.any():
            break
    if "compound_key" in df.columns:
        df["compound_key"] = df["compound_key"].astype("string")

    log_value: pd.Series | None = None
    if "log_value" in df.columns:
        log_value = pd.Series(
            pd.to_numeric(df["log_value"], errors="coerce"),
            index=df.index,
            dtype="Float64",
        )
    elif "pchembl_value" in df.columns:
        log_value = pd.Series(
            pd.to_numeric(df["pchembl_value"], errors="coerce"),
            index=df.index,
            dtype="Float64",
        )
        filled.add("log_value")

    if log_value is not None:
        df["log_value"] = log_value

    if (
        log_value is not None
        and "standard_value" in df.columns
        and "standard_units" in df.columns
    ):
        std_value = pd.to_numeric(df["standard_value"], errors="coerce")
        units = df["standard_units"].astype("string").str.strip().str.lower()
        factors = units.map({
            "m": 1.0,
            "mm": 1e-3,
            "um": 1e-6,
            "µm": 1e-6,
            "nm": 1e-9,
            "pm": 1e-12,
        })
        factors = pd.to_numeric(factors, errors="coerce")
        molar = std_value * factors
        computed = pd.Series(pd.NA, index=df.index, dtype="Float64")
        valid = molar.notna() & (molar > 0)
        if valid.any():
            computed_values = -np.log10(molar[valid].astype("float64"))
            computed.loc[valid] = pd.Series(computed_values, index=df.index[valid]).astype(
                "Float64"
            )
        mask = log_value.isna() & computed.notna()
        if mask.any():
            log_value = log_value.mask(mask, computed)
            df["log_value"] = log_value

    bool_defaults = {
        "multmol_assay",
        "approx_cited_activity",
        "shuffled_cit",
        "exact_cited_activity",
        "higly_correlated_cit",
        "review_doc",
        "rounded_data_citation",
    }
    for column in bool_defaults:
        if column not in df.columns:
            df[column] = pd.Series(False, index=df.index, dtype="boolean")
            filled.add(column)

    float_defaults = {"original_activity_approx", "original_activity_exact"}
    for column in float_defaults:
        if column not in df.columns:
            df[column] = pd.Series(pd.NA, index=df.index, dtype="Float64")
            filled.add(column)

    if "salt_chembl_id" in df.columns and "molecule_chembl_id" in df.columns:
        salt_series = df["salt_chembl_id"].astype("string")
        molecule_series = df["molecule_chembl_id"].astype("string")
        parent_mask = (
            df["parent_molecule_chembl_id"].notna()
            if "parent_molecule_chembl_id" in df.columns
            else pd.Series(False, index=df.index)
        )
        same = salt_series.eq(molecule_series) & parent_mask
        same = same.fillna(False)
        if same.any():
            df.loc[same, "salt_chembl_id"] = pd.NA
            filled.add("salt_chembl_id")

    if "salt_chembl_id" not in df.columns:
        df["salt_chembl_id"] = pd.Series(pd.NA, index=df.index, dtype="string")
        filled.add("salt_chembl_id")
    else:
        df["salt_chembl_id"] = df["salt_chembl_id"].astype("string")

    if "nstereo" not in df.columns:
        df["nstereo"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
        filled.add("nstereo")

    return df, filled


def _format_log_value_for_export(value: object) -> object:
    if pd.isna(value):  # type: ignore[arg-type]
        return pd.NA
    text = format(float(value), ".15g")
    if "e" not in text.lower() and "." not in text:
        text = f"{text}.0"
    return text


def _transform_activity_frame(
    frame: pd.DataFrame,
    *,
    dictionary_root: Path,
    targets_override: Path | None,
) -> pd.DataFrame:
    df, ensured = _ensure_required_input_columns(frame)
    df, augmented = _augment_activity_frame(df)
    filled = ensured | augmented
    missing = _REQUIRED_COLUMNS - set(df.columns)
    if missing:
        available = ", ".join(sorted(df.columns))
        missing_list = ", ".join(sorted(missing))
        raise ActivityExtendedError(
            "Activity export missing required columns: "
            f"{missing_list}. Available columns: {available}"
        )

    working, filled = _augment_activity_frame(df)

    if filled:
        logger.warning(
            "activity_extended_missing_columns_filled",
            columns=sorted(filled),
        )

    df = _prepare_unknown_chirality(df)
    df = _apply_multimol_logic(df)
    df = _rename_columns(df)
    df = _drop_unused_columns(df)
    df = _compute_citation_flags(df)
    df = _annotate_high_citation(df, dictionary_root)
    df = _merge_target_metadata(df, dictionary_root=dictionary_root, targets_override=targets_override)
    df = _select_and_cast(df)
    return df


def _derive_output_path(input_path: Path) -> Path:
    base_name = helpers.normalise_export_basename(input_path)
    candidate = base_name
    while candidate.startswith("extended."):
        candidate = candidate[len("extended.") :]

    match = _FILENAME_RE.match(candidate)
    if match is not None:
        stamp = match.group(1)
        final_name = f"extended.output.activity_{stamp}.csv"
    else:
        logger.warning("activity_extended_unexpected_basename", basename=base_name)
        final_name = f"extended.{candidate}"

    return input_path.with_name(final_name)


def process_activity_extended(
    search_dir: Path | str | None = None,
    *,
    input_path: Path | str | None = None,
    dictionary_dir: Path | str | None = None,
    targets_csv: Path | str | None = None,
    base_dir: Path | str | None = None,
) -> Path:
    """Create the extended activity export mirroring the Power Query workflow.

    Parameters
    ----------
    search_dir:
        Directory containing ``output.activity_<stamp>.csv`` exports.  When
        omitted, ``data/output`` is scanned unless ``input_path`` is provided.
    input_path:
        Explicit path to an activity export. When supplied the search directory
        scan is skipped and ``input_path`` is used directly. This is useful when
        the caller already knows the exact export path (for example when the
        file was just produced by the pipeline) or when the export follows a
        non-standard naming convention.
    dictionary_dir:
        Root directory with bundled dictionary CSVs.  Defaults to
        ``config/dictionary``.
    targets_csv:
        Optional explicit path to ``targets_type.csv`` overriding the dictionary
        lookup.

    base_dir:
        Deprecated alias for ``search_dir`` kept for backwards compatibility.

    Returns
    -------
    pathlib.Path
        Path to the generated ``extended.output.activity_<stamp>.csv`` file.
    """

    if base_dir is not None:
        if search_dir is not None:
            raise TypeError(
                "process_activity_extended() received both 'search_dir' and 'base_dir'. "
                "Use 'search_dir' only."
            )
        warnings.warn(
            "'base_dir' is deprecated; pass 'search_dir' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        search_dir = base_dir

    explicit_input: Path | None = Path(input_path) if input_path is not None else None
    if explicit_input is not None:
        if explicit_input.is_dir():
            raise ActivityExtendedError(
                f"input_path must point to a file, received directory: {explicit_input!s}"
            )
        if not explicit_input.exists():
            raise ActivityExtendedError(f"Activity export not found: {explicit_input!s}")

    resolved_search_dir = (
        Path(search_dir)
        if search_dir is not None
        else (explicit_input.parent if explicit_input is not None else _current_default_search_dir())
    )
    if not resolved_search_dir.exists():
        raise ActivityExtendedError(f"Search directory does not exist: {resolved_search_dir!s}")
    if not resolved_search_dir.is_dir():
        raise ActivityExtendedError(
            f"Search directory is not a directory: {resolved_search_dir!s}"
        )

    input_path = explicit_input or _latest_activity_export(resolved_search_dir)
    dictionary_root = _resolve_dictionary_root(Path(dictionary_dir) if dictionary_dir is not None else None)

    frame = helpers.read_csv_with_fallbacks(input_path)
    processed = _transform_activity_frame(
        frame,
        dictionary_root=dictionary_root,
        targets_override=Path(targets_csv) if targets_csv is not None else None,
    )

    output_path = _derive_output_path(input_path)
    logger.info(
        "activity_extended_dataframe_shape",
        rows=processed.shape[0],
        columns=processed.shape[1],
    )
    logger.info(
        "activity_extended_columns",
        columns=sorted(processed.columns.tolist()),
    )
    for column in ("activity_chembl_id", "saltform_id", "pA_value"):
        non_null = int(processed[column].notna().sum()) if column in processed.columns else 0
        logger.info(
            "activity_extended_non_null_counts",
            column=column,
            non_null=non_null,
            total=len(processed),
        )
    logger.info("activity_extended_saving", path=str(output_path))
    output_frame = processed.copy()
    if "pA_value" in output_frame.columns:
        output_frame["pA_value"] = (
            output_frame["pA_value"].map(_format_log_value_for_export).astype("string")
        )
    helpers.write_csv(output_frame, output_path, columns=_FINAL_COLUMN_ORDER)
    logger.info(
        "activity_extended_saved",
        path=str(output_path),
        rows=len(processed),
    )
    return output_path


__all__ = ["process_activity_extended", "ActivityExtendedError"]

