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

import functools
import json
import re
import sys
import warnings
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from config.paths import DICTIONARY_DIR
from library.common.log import logger

from . import helpers

_DEFAULT_SEARCH_DIR = Path("data/output")
_DEFAULT_DICTIONARY_DIR = DICTIONARY_DIR
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

_OPTIONAL_EMPTY_BACKFILL_COLUMNS: frozenset[str] = frozenset(
    {
        "multmol_assay",
        "original_activity_approx",
        "original_activity_exact",
    }
)

_NA_MARKERS: tuple[str, ...] = ("[#N/A]",)

_ACTIVITY_INPUT_SCHEMA: Mapping[str, str] = {
    column: (
        "Logical"
        if dtype == "boolean"
        else (
            "Int64" if dtype == "Int64" else "Float64" if dtype == "Float64" else "Text"
        )
    )
    for column, dtype in _REQUIRED_COLUMN_DTYPES.items()
    if dtype in {"string", "boolean", "Int64", "Float64"}
}


def _normalise_column_key(name: str) -> str:
    """Return a canonical key used to match column aliases."""

    return re.sub(r"[^0-9a-z]+", "_", str(name).strip().lower()).strip("_")


def _lookup_column_name(frame: pd.DataFrame, *candidates: str) -> str | None:
    """Return the actual column name matching one of ``candidates``."""

    normalised = {_normalise_column_key(column): column for column in frame.columns}
    for candidate in candidates:
        key = _normalise_column_key(candidate)
        if key in normalised:
            return normalised[key]
    return None


_TARGET_METADATA_READ_SCHEMA: Mapping[str, str] = {
    "target_chembl_id": "Text",
    "target_sort_order": "Text",
    "multifunctional_enzyme": "Logical",
    "IUPHAR_class": "Text",
    "IUPHAR_subclass": "Text",
    "genus": "Text",
    "superkingdom": "Text",
    "phylum": "Text",
    "taxon_id": "Int64",
    "gene_index": "Text",
    "taxon_index": "Text",
    "unicellular_organism": "Logical",
}

_TARGET_METADATA_SCHEMA: Mapping[str, str] = {
    "target_chembl_id": "Text",
    "sortorder.target": "Text",
    "multifunctional_enzyme": "Logical",
    "IUPHAR_class": "Text",
    "IUPHAR_subclass": "Text",
    "genus": "Text",
    "superkingdom": "Text",
    "phylum": "Text",
    "taxon_id": "Int64",
    "gene_index": "Text",
    "taxon_index": "Text",
    "unicellular_organism": "Logical",
}

_CITATION_FRACTION_SCHEMA: Mapping[str, str] = {
    "N": "Int64",
    "K_min_significant": "Int64",
}

_REQUIRED_COLUMN_FALLBACKS: Mapping[str, Callable[[pd.DataFrame], pd.Series | None]] = {
    "activity_chembl_id": lambda frame: frame.get("activity_id"),
    "salt_chembl_id": (
        lambda frame: (
            frame.get("molecule_chembl_id")
            if "parent_molecule_chembl_id" not in frame.columns
            else None
        )
    ),
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

_DEDUPE_SORT_ORDER: tuple[str, ...] = (
    "completed",
    "activity_id",
    "activity_chembl_id",
    "assay_chembl_id",
    "document_chembl_id",
    "standard_type",
    "standard_value",
)

_DEDUPE_SUBSET_KEYS: tuple[str, ...] = (
    "activity_id",
    "assay_chembl_id",
    "document_chembl_id",
    "standard_type",
    "standard_value",
)

_FINAL_COLUMN_ORDER: tuple[str, ...] = (
    "activity_chembl_id",
    "saltform_id",
    "molecule_chembl_id",
    "molecule_chembl_id.1",
    "target_chembl_id",
    "assay_chembl_id",
    "document_chembl_id",
    "completed",
    "bao_endpoint",
    "standard_type",
    "standard_value",
    "pA_value",
    "bao_format",
    "compound_key",
    "compound_name",
    "standard_inchi_skeleton",
    "unknown_chirality",
    "multmol_assay",
    "assay_with_same_target",
    "exact_data_citation",
    "higly_correlated_assay",
    "shuffled_assay",
    "review",
    "rounded_data_citation",
    "allosteric",
    "nam",
    "pam",
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
    "sortorder.target",
)

_TARGET_COLUMNS: Sequence[str] = (
    "target_chembl_id",
    "sortorder.target",
    "multifunctional_enzyme",
    "IUPHAR_class",
    "IUPHAR_subclass",
    "genus",
    "superkingdom",
    "phylum",
    "taxon_id",
    "gene_index",
    "taxon_index",
    "unicellular_organism",
)

_UNICELLULAR_NORMALISATIONS: Mapping[str, bool] = {
    "unicellular": True,
    "unicellular organism": True,
    "multicellular": False,
    "multicellular organism": False,
}

_NA_TYPE = type(pd.NA)


class ActivityExtendedError(RuntimeError):
    """Raised when the activity extended post-processing cannot proceed."""


def _current_default_search_dir() -> Path:
    package = sys.modules.get(__name__)
    if package is not None and hasattr(package, "_DEFAULT_SEARCH_DIR"):
        override = package._DEFAULT_SEARCH_DIR
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


def _resolve_dictionary_root(dictionary_dir: Path | str | None) -> Path:
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
    return helpers.read_csv_strict(
        candidate,
        encoding=helpers.ENCODING_FALLBACKS,
        dtype_map=_CITATION_FRACTION_SCHEMA,
        na_values=_NA_MARKERS,
    )


def _resolve_targets_path(dictionary_root: Path, override: Path | None) -> Path:
    if override is not None:
        path = Path(override)
        if not path.exists():
            raise ActivityExtendedError(
                f"targets_type.csv override not found: {path!s}"
            )
        return path

    candidates = [
        dictionary_root / "targets_type.csv",
        dictionary_root / "_target" / "targets_type.csv",
        dictionary_root / "_Target" / "targets_type.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    formatted = (
        " or ".join(f"'{path}'" for path in candidates[:-1]) + f" or '{candidates[-1]}'"
    )
    raise ActivityExtendedError(
        "targets_type.csv not found in the provided dictionary directory. Expected at "
        + formatted
    )


def _load_target_metadata(path: Path) -> pd.DataFrame:
    frame = helpers.read_csv_strict(
        path,
        encoding="cp1252",
        dtype_map=_TARGET_METADATA_READ_SCHEMA,
        na_values=_NA_MARKERS,
    )
    alias = _lookup_column_name(frame, "unicellular_organism", "unicellular organism")
    if alias is not None and alias != "unicellular_organism":
        frame = frame.rename(columns={alias: "unicellular_organism"})
    if "unicellular_organism" not in frame.columns:
        source_column = _lookup_column_name(
            frame, "type", "organism_type", "organism type"
        )
        if source_column is not None:
            source = frame[source_column].astype("string")
            normalised = source.str.strip().str.lower()
            inferred = normalised == "unicellular organism"
            frame["unicellular_organism"] = inferred.astype("boolean")
        else:
            frame["unicellular_organism"] = pd.Series(
                pd.NA, index=frame.index, dtype="boolean"
            )
    frame = frame.rename(columns={"target_sort_order": "sortorder.target"})

    if "unicellular_organism" not in frame.columns and "organism_type" in frame.columns:
        frame = frame.copy()

        unmapped_values: set[str] = set()

        def _normalise_organism(value: object) -> bool | _NA_TYPE:
            if pd.isna(value):
                return pd.NA
            text = str(value).strip()
            if not text:
                return pd.NA
            resolved = _UNICELLULAR_NORMALISATIONS.get(text.casefold())
            if resolved is None:
                unmapped_values.add(text)
                return pd.NA
            return resolved

        inferred = frame["organism_type"].map(_normalise_organism).astype("boolean")
        if unmapped_values:
            logger.warning(
                "targets_type.csv contains organism_type values without boolean mapping: %s",
                sorted(unmapped_values),
            )
        frame["unicellular_organism"] = inferred

    missing = [column for column in _TARGET_COLUMNS if column not in frame.columns]
    if missing:
        raise ActivityExtendedError(
            "targets_type.csv missing expected columns: " + ", ".join(sorted(missing))
        )
    return frame.loc[:, list(_TARGET_COLUMNS)]


def _load_document_lookup(dictionary_root: Path) -> pd.DataFrame:
    candidate = dictionary_root / "_document" / "document.csv"
    if not candidate.exists():
        raise ActivityExtendedError(
            "document.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
        )

    frame: pd.DataFrame | None = None
    missing_errors: list[str] = []
    for sep in helpers.CSV_SEPARATORS:
        try:
            candidate_frame = helpers.read_csv_with_fallbacks(candidate, sep=sep)
        except Exception as exc:
            missing_errors.append(f"sep={sep!r}: {exc!s}")
            continue

        resolved_columns: dict[str, str] = {}
        for column in candidate_frame.columns:
            normalised = column.strip().lower()
            if normalised == "document_chembl_id":
                resolved_columns[column] = "document_chembl_id"
            elif normalised in {
                "chembl.document_chembl_id",
                "chembl.document chembl id",
            }:
                resolved_columns[column] = "document_chembl_id"

        if resolved_columns:
            candidate_frame = candidate_frame.rename(columns=resolved_columns)
            alias_source = [
                src
                for src, dst in resolved_columns.items()
                if dst == "document_chembl_id" and src != dst
            ]
            if alias_source:
                logger.warning(
                    "Renamed document identifier column(s) %s to 'document_chembl_id'. "
                    "Update the dictionary export to align with the expected schema.",
                    alias_source,
                )

        if "document_chembl_id" in candidate_frame.columns:
            frame = candidate_frame
            break

        missing_errors.append(
            "sep={}: available columns -> {}".format(
                repr(sep), ", ".join(candidate_frame.columns) or "<none>"
            )
        )

    if frame is None:
        detail = "; ".join(missing_errors)
        raise ActivityExtendedError(
            "document.csv missing required 'document_chembl_id' column; unable to join document metadata"
            + (f". {detail}" if detail else "")
        )

    columns = ["document_chembl_id", "completed", "review"]
    for column in columns:
        if column not in frame.columns:
            frame[column] = pd.Series(pd.NA, dtype="string")
    result = frame.loc[:, columns].drop_duplicates(subset=["document_chembl_id"])
    return result


def _load_assay_lookup(dictionary_root: Path) -> pd.DataFrame:
    candidate = dictionary_root / "_assay" / "assay.csv"
    if not candidate.exists():
        raise ActivityExtendedError(
            "assay.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
        )
    frame: pd.DataFrame | None = None
    for sep in (",", "\t"):
        try:
            candidate_frame = helpers.read_csv_with_fallbacks(candidate, sep=sep)
        except Exception:
            continue
        if "assay_chembl_id" in candidate_frame.columns:
            frame = candidate_frame
            break
    if frame is None:
        raise ActivityExtendedError(
            "assay.csv missing required 'assay_chembl_id' column; unable to join assay metadata"
        )
    if "assay_with_same_target" not in frame.columns:
        frame["assay_with_same_target"] = pd.Series(pd.NA, dtype="string")
    result = frame.loc[:, ["assay_chembl_id", "assay_with_same_target"]]
    return result.drop_duplicates(subset=["assay_chembl_id"])


def _load_testitem_lookup(dictionary_root: Path) -> pd.DataFrame:
    candidate = dictionary_root / "_testitem" / "testitem.csv"
    if not candidate.exists():
        raise ActivityExtendedError(
            "testitem.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
        )
    frame: pd.DataFrame | None = None
    for sep in (",", "\t"):
        try:
            candidate_frame = helpers.read_csv_with_fallbacks(candidate, sep=sep)
        except Exception:
            continue
        if "molecule_chembl_id" in candidate_frame.columns:
            frame = candidate_frame
            break
    if frame is None:
        raise ActivityExtendedError(
            "testitem.csv missing required 'molecule_chembl_id' column; unable to join testitem metadata"
        )
    if "standard_inchi_skeleton" not in frame.columns:
        frame["standard_inchi_skeleton"] = pd.Series(pd.NA, dtype="string")
    result = frame.loc[:, ["molecule_chembl_id", "standard_inchi_skeleton"]]
    return result.drop_duplicates(subset=["molecule_chembl_id"])


@functools.cache
def _load_parent_lookup_cached(dictionary_root: str) -> pd.Series:
    root_path = Path(dictionary_root)
    candidate = root_path / "_testitem" / "molecule_hierarchy.csv"
    if not candidate.exists():
        raise ActivityExtendedError(
            "molecule_hierarchy.csv not found; expected at "
            f"'{candidate}'. Provide dictionary_dir pointing to the bundled dictionaries."
        )
    try:
        frame = helpers.read_csv_with_fallbacks(candidate)
    except Exception as exc:  # pragma: no cover - defensive
        raise ActivityExtendedError(
            "unable to read molecule_hierarchy.csv; verify the CSV format"
        ) from exc

    required = {"molecule_chembl_id", "parent_molecule_chembl_id"}
    missing = required - set(frame.columns)
    if missing:
        missing_list = ", ".join(sorted(missing))
        raise ActivityExtendedError(
            "molecule_hierarchy.csv missing required columns: " + missing_list
        )

    subset = frame.loc[:, ["molecule_chembl_id", "parent_molecule_chembl_id"]].copy()
    subset["molecule_chembl_id"] = subset["molecule_chembl_id"].astype("string")
    subset["parent_molecule_chembl_id"] = subset["parent_molecule_chembl_id"].astype(
        "string"
    )
    subset = subset.dropna(subset=["molecule_chembl_id"])
    subset = subset.loc[subset["molecule_chembl_id"].str.strip().ne("")]
    subset = subset.drop_duplicates(subset=["molecule_chembl_id"], keep="first")
    return subset.set_index("molecule_chembl_id")["parent_molecule_chembl_id"]


def _load_parent_lookup(dictionary_root: Path) -> pd.Series:
    resolved_root = dictionary_root.resolve()
    return _load_parent_lookup_cached(str(resolved_root))


def _insert_columns_after(
    df: pd.DataFrame, anchor: str, columns: Sequence[str]
) -> pd.DataFrame:
    present = [column for column in columns if column in df.columns]
    if not present:
        return df
    base_order = [column for column in df.columns if column not in present]
    try:
        index = base_order.index(anchor) + 1
    except ValueError:
        index = len(base_order)
    new_order = base_order[:index] + present + base_order[index:]
    return df.loc[:, new_order]


def _log_join_statistics(event: str, indicator: pd.Series) -> None:
    matched = int((indicator == "both").sum())
    missing = int((indicator == "left_only").sum())
    logger.info(event, matched=matched, missing=missing)


def _safe_to_bool(series: pd.Series, column: str) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise ActivityExtendedError(
            f"column '{column}' has duplicate entries; expected a Series"
        )

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
        raise ActivityExtendedError(
            f"column '{column}' has duplicate entries; expected a Series"
        )
    try:
        return pd.to_numeric(series, errors="raise").astype("Int64")
    except Exception as exc:  # pragma: no cover - defensive downgrade
        logger.warning("dtype_int_conversion_failed", column=column, error=str(exc))
        return series.astype("string")


def _prepare_unknown_chirality(frame: pd.DataFrame) -> pd.DataFrame:
    df = frame.copy()
    if "nstereo" in df.columns:
        df["unknown_chirality"] = (
            _safe_to_int(df["nstereo"], "nstereo").ne(1).fillna(True)
        )
        df.drop(columns=["nstereo"], inplace=True)
    else:
        df["unknown_chirality"] = pd.Series(True, index=df.index, dtype="boolean")
    return df


def _apply_multimol_logic(df: pd.DataFrame) -> pd.DataFrame:
    missing = set(_GROUP_KEY_COLUMNS) - set(df.columns)
    if missing:
        raise ActivityExtendedError(
            "activity table missing columns for multimol grouping: "
            + ", ".join(sorted(missing))
        )
    df = helpers.sort_power_query(df, _GROUP_KEY_COLUMNS)
    counts = (
        df.groupby(list(_GROUP_KEY_COLUMNS), dropna=False)
        .size()
        .rename("Count")
        .reset_index()
    )
    counts = helpers.sort_power_query(counts, _GROUP_KEY_COLUMNS)
    merged = df.merge(counts, on=list(_GROUP_KEY_COLUMNS), how="left")
    mask = (
        merged["unknown_chirality"].fillna(True).eq(False)
        & merged["multmol_assay"].isna()
        & (merged["Count"] > 1)
    )
    duplicated_assays = set(merged.loc[mask, "assay_chembl_id"].dropna().astype(str))
    merged["multimol_assay_same"] = merged["assay_chembl_id"].isin(duplicated_assays)

    multmol_series = _safe_to_bool(merged["multmol_assay"], "multmol_assay").fillna(
        False
    )
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


def _ensure_required_input_columns(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str]]:
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


def _ensure_compound_key_sources(
    frame: pd.DataFrame, *, dictionary_root: Path
) -> tuple[pd.DataFrame, set[str]]:
    df = frame.copy()
    filled: set[str] = set()

    if "molecule_chembl_id" not in df.columns:
        df["molecule_chembl_id"] = pd.Series(pd.NA, index=df.index, dtype="string")
        filled.add("molecule_chembl_id")

    if "parent_molecule_chembl_id" not in df.columns:
        df["parent_molecule_chembl_id"] = pd.Series(
            pd.NA, index=df.index, dtype="string"
        )
        filled.add("parent_molecule_chembl_id")

    if df.empty:
        return df, filled

    parent_missing = _string_missing_mask(df["parent_molecule_chembl_id"])
    molecule_available = ~_string_missing_mask(df["molecule_chembl_id"])
    needs_parent = parent_missing & molecule_available
    if needs_parent.any():
        lookup = _load_parent_lookup(dictionary_root)
        resolved = df.loc[needs_parent, "molecule_chembl_id"].map(lookup)
        resolved = resolved.dropna()
        if not resolved.empty:
            df.loc[resolved.index, "parent_molecule_chembl_id"] = resolved
            filled.add("parent_molecule_chembl_id")

    parent_missing = _string_missing_mask(df["parent_molecule_chembl_id"])
    molecule_missing = _string_missing_mask(df["molecule_chembl_id"])
    unresolved = parent_missing & molecule_missing
    if unresolved.any():
        raise ActivityExtendedError(
            "Unable to derive compound_key: missing both molecule_chembl_id and "
            "parent_molecule_chembl_id for one or more rows"
        )

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


def _normalise_activity_properties_text(value: object) -> str | None:
    if value is None:
        return None
    if isinstance(value, float | np.floating) and np.isnan(value):
        return None
    if value is pd.NA:
        return None
    text = str(value) if not isinstance(value, str) else value
    text = text.strip()
    if not text:
        return None
    if text.lower() in {"nan", "none", "null"}:
        return None
    if text.startswith('"') and text.endswith('"'):
        inner = text[1:-1]
        if inner.startswith("{") and inner.endswith("}"):
            text = inner
    return text


def _load_activity_properties_json(text: str) -> Mapping[str, object] | None:
    candidates = [text]
    if '""' in text:
        candidates.append(text.replace('""', '"'))
    if text.startswith('"') and text.endswith('"'):
        inner = text[1:-1]
        candidates.append(inner)
        candidates.append(inner.replace('""', '"'))

    expanded: list[str] = []
    for candidate in candidates:
        replaced = (
            candidate.replace("True", "true")
            .replace("False", "false")
            .replace("None", "null")
        )
        expanded.extend([candidate, replaced, replaced.replace("'", '"')])

    seen: set[str] = set()
    for candidate in expanded:
        stripped = candidate.strip()
        if not stripped or stripped in seen:
            continue
        seen.add(stripped)
        try:
            loaded = json.loads(stripped)
        except json.JSONDecodeError:
            continue
        if isinstance(loaded, Mapping):
            return loaded
    return None


def _coerce_activity_property_flag(value: object) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, int | np.integer):
        if value in (0, 1):
            return bool(value)
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes", "y"}:
            return True
        if lowered in {"false", "0", "no", "n"}:
            return False
    return None


def _extract_activity_properties_flags(df: pd.DataFrame) -> pd.DataFrame:
    result = df.copy()
    index = result.index
    defaults = {
        column: pd.Series(pd.NA, index=index, dtype="boolean")
        for column in ("allosteric", "nam", "pam")
    }
    flags = pd.DataFrame(defaults)
    series = result.get("activity_properties")
    if series is None:
        for column in flags.columns:
            result[column] = flags[column]
        return result

    mapping = {
        "allosteric": ("allosteric",),
        "nam": ("negative", "nam"),
        "pam": ("positive", "pam"),
    }

    for idx, raw in series.items():
        text = _normalise_activity_properties_text(raw)
        if not text:
            continue
        payload = _load_activity_properties_json(text)
        if payload is None:
            continue
        feature_source: Mapping[str, object] | None = None
        features = (
            payload.get("effect_features") if isinstance(payload, Mapping) else None
        )
        if isinstance(features, Mapping):
            feature_source = features
        elif isinstance(payload, Mapping):
            feature_source = payload
        else:
            feature_source = {}

        for column, keys in mapping.items():
            if feature_source is None:
                continue
            value: object | None = None
            for key in keys:
                if key in feature_source:
                    value = feature_source[key]
                    break
            coerced = _coerce_activity_property_flag(value)
            flags.at[idx, column] = False if coerced is None else coerced

    for column in flags.columns:
        result[column] = flags[column]
    return result


def _compute_citation_flags(df: pd.DataFrame) -> pd.DataFrame:
    bool_columns = [
        "exact_data_citation",
        "higly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
    ]
    prepared: dict[str, pd.Series] = {}
    for column in bool_columns:
        if column in df.columns:
            prepared[column] = _safe_to_bool(df[column], column).fillna(False)
        else:
            prepared[column] = pd.Series(False, index=df.index, dtype="boolean")
    mask = pd.DataFrame(prepared).any(axis=1).astype("boolean")
    result = df.copy()
    result["is_citation"] = mask
    return result


def _annotate_high_citation(df: pd.DataFrame, dictionary_root: Path) -> pd.DataFrame:
    converted = helpers.sort_power_query(
        df, ("document_chembl_id", "activity_chembl_id")
    )
    counts = (
        converted.groupby("document_chembl_id")["is_citation"]
        .agg(n_citation="sum", n_non_citation=lambda s: (~s).sum())
        .reset_index()
    )
    counts["N"] = counts["n_citation"] + counts["n_non_citation"]
    counts = counts[(counts["n_citation"] > 0) & (counts["n_non_citation"] > 0)]

    citation_fraction = _load_citation_fraction(dictionary_root)
    counts = helpers.sort_power_query(counts, ("document_chembl_id",))
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


def _merge_document_metadata(df: pd.DataFrame, dictionary_root: Path) -> pd.DataFrame:
    lookup = _load_document_lookup(dictionary_root)
    merged = df.merge(
        lookup,
        on="document_chembl_id",
        how="left",
        indicator="_merge_document",
    )
    _log_join_statistics("activity_extended_document_join", merged["_merge_document"])
    merged.drop(columns=["_merge_document"], inplace=True)

    if "completed" in merged.columns:
        merged["completed"] = merged["completed"].astype("string")
    else:
        merged["completed"] = pd.Series(pd.NA, index=merged.index, dtype="string")

    if "review" in merged.columns:
        merged["review"] = _safe_to_bool(merged["review"], "review")
    else:
        merged["review"] = pd.Series(pd.NA, index=merged.index, dtype="boolean")

    merged = _insert_columns_after(
        merged, "document_chembl_id", ("completed", "review")
    )
    return merged


def _merge_assay_metadata(df: pd.DataFrame, dictionary_root: Path) -> pd.DataFrame:
    lookup = _load_assay_lookup(dictionary_root)
    merged = df.merge(
        lookup,
        on="assay_chembl_id",
        how="left",
        indicator="_merge_assay",
        suffixes=("", "_assay"),
    )
    _log_join_statistics("activity_extended_assay_join", merged["_merge_assay"])
    merged.drop(columns=["_merge_assay"], inplace=True)

    if "assay_with_same_target" in merged.columns:
        merged["assay_with_same_target"] = _safe_to_int(
            merged["assay_with_same_target"], "assay_with_same_target"
        )
    else:
        merged["assay_with_same_target"] = pd.Series(
            pd.NA, index=merged.index, dtype="Int64"
        )

    merged = _insert_columns_after(merged, "multmol_assay", ("assay_with_same_target",))
    return merged


def _merge_testitem_metadata(df: pd.DataFrame, dictionary_root: Path) -> pd.DataFrame:
    lookup = _load_testitem_lookup(dictionary_root)
    lookup = lookup.assign(_molecule_chembl_id_expanded=lookup["molecule_chembl_id"])
    merged = df.merge(
        lookup,
        on="molecule_chembl_id",
        how="left",
        indicator="_merge_testitem",
        suffixes=("", "_testitem"),
    )
    _log_join_statistics("activity_extended_testitem_join", merged["_merge_testitem"])
    merged.drop(columns=["_merge_testitem"], inplace=True)

    right_id_column = "_molecule_chembl_id_expanded"
    if right_id_column in merged.columns:
        new_ids = merged[right_id_column].astype("string")
        merged.drop(columns=[right_id_column], inplace=True)
    else:
        new_ids = pd.Series(pd.NA, index=merged.index, dtype="string")
    merged["molecule_chembl_id.1"] = new_ids

    if "standard_inchi_skeleton" in merged.columns:
        left_inchi = merged["standard_inchi_skeleton"].astype("string")
        merged.drop(columns=["standard_inchi_skeleton"], inplace=True)
    else:
        left_inchi = pd.Series(pd.NA, index=merged.index, dtype="string")

    right_inchi_column = "standard_inchi_skeleton_testitem"
    if right_inchi_column in merged.columns:
        right_inchi = merged[right_inchi_column].astype("string")
        merged.drop(columns=[right_inchi_column], inplace=True)
    else:
        right_inchi = pd.Series(pd.NA, index=merged.index, dtype="string")
    merged["standard_inchi_skeleton"] = right_inchi.fillna(left_inchi)

    merged["molecule_chembl_id.1"] = merged["molecule_chembl_id.1"].astype("string")
    merged["standard_inchi_skeleton"] = merged["standard_inchi_skeleton"].astype(
        "string"
    )

    merged = _insert_columns_after(
        merged, "molecule_chembl_id", ("molecule_chembl_id.1",)
    )
    merged = _insert_columns_after(
        merged, "compound_name", ("standard_inchi_skeleton",)
    )
    return merged


def _merge_target_metadata(
    df: pd.DataFrame,
    *,
    dictionary_root: Path,
    targets_override: Path | None,
) -> pd.DataFrame:
    targets_path = _resolve_targets_path(dictionary_root, targets_override)
    targets = _load_target_metadata(targets_path)
    df = helpers.sort_power_query(df, ("target_chembl_id", "assay_chembl_id"))
    targets = helpers.sort_power_query(targets, ("target_chembl_id",))
    merged = df.merge(targets, on="target_chembl_id", how="left")
    merged = merged.loc[:, ~merged.columns.duplicated()]

    merged = helpers.coerce_types(merged, _TARGET_METADATA_SCHEMA)

    for column in ("multifunctional_enzyme", "unicellular_organism"):
        if column in merged.columns:
            merged[column] = _safe_to_bool(merged[column], column)

    reorder_plan: tuple[tuple[str, tuple[str, ...]], ...] = (
        ("original_activity_exact", ("IUPHAR_class", "IUPHAR_subclass")),
        ("IUPHAR_subclass", ("unicellular_organism", "multifunctional_enzyme")),
        ("multifunctional_enzyme", ("gene_index", "taxon_index")),
        ("taxon_index", ("sortorder.target",)),
    )
    for anchor, columns in reorder_plan:
        merged = _insert_columns_after(merged, anchor, columns)

    drop_columns = [
        column
        for column in ("genus", "superkingdom", "phylum", "taxon_id")
        if column in merged.columns
    ]
    if drop_columns:
        merged.drop(columns=drop_columns, inplace=True)

    return merged


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
        "allosteric",
        "nam",
        "pam",
        "high_citation_rate",
        "is_citation",
        "unicellular_organism",
        "multifunctional_enzyme",
    ]
    for column in bool_columns:
        if column in result.columns:
            result[column] = _safe_to_bool(result[column], column)

    def _format_numeric(value: object) -> object:
        if pd.isna(value):
            return pd.NA
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            return value
        if float(numeric).is_integer():
            return f"{numeric:.1f}"
        text = str(numeric)
        return text.rstrip("0").rstrip(".") if "." in text else text

    for column in ("standard_value", "pA_value"):
        if column in result.columns:
            result[column] = result[column].map(_format_numeric).astype("string")
    return result


def dedupe_final(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` sorted and deduplicated using Power Query key rules."""

    if df.empty:
        logger.info("activity_extended_deduplicated", removed=0, remaining=0)
        return df.copy()

    sort_candidates: list[str] = []
    for column in _DEDUPE_SORT_ORDER:
        if column == "activity_id":
            options = ("activity_id", "activity_chembl_id")
        else:
            options = (column,)
        for option in options:
            if option in df.columns and option not in sort_candidates:
                sort_candidates.append(option)
                break

    sorted_df = (
        helpers.sort_power_query(df, sort_candidates) if sort_candidates else df.copy()
    )

    subset: list[str] = []
    missing: list[str] = []
    for column in _DEDUPE_SUBSET_KEYS:
        options = (
            ("activity_id", "activity_chembl_id")
            if column == "activity_id"
            else (column,)
        )
        selected = next(
            (candidate for candidate in options if candidate in sorted_df.columns), None
        )
        if selected is None:
            missing.append(column)
        else:
            if selected not in subset:
                subset.append(selected)

    if missing:
        raise ActivityExtendedError(
            "activity table missing columns required for deduplication: "
            + ", ".join(sorted(missing))
        )

    deduped = sorted_df.drop_duplicates(subset=subset, keep="first").reset_index(
        drop=True
    )
    removed = int(len(sorted_df) - len(deduped))
    logger.info(
        "activity_extended_deduplicated", removed=removed, remaining=len(deduped)
    )
    return deduped


def _string_missing_mask(series: pd.Series) -> pd.Series:
    values = series.astype("string")
    stripped = values.str.strip()
    return values.isna() | stripped.fillna("").eq("")


def _augment_activity_frame(frame: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    df = frame.copy()
    filled: set[str] = set()

    if "activity_chembl_id" not in df.columns:
        if "activity_id" in df.columns:
            df["activity_chembl_id"] = df["activity_id"].astype("string")
            filled.add("activity_chembl_id")
    else:
        mask = _string_missing_mask(df["activity_chembl_id"])
        if mask.any() and "activity_id" in df.columns:
            df.loc[mask, "activity_chembl_id"] = df.loc[mask, "activity_id"].astype(
                "string"
            )
            filled.add("activity_chembl_id")

    if "compound_name" not in df.columns:
        if "molecule_pref_name" in df.columns:
            df["compound_name"] = df["molecule_pref_name"].astype("string")
            filled.add("compound_name")
    else:
        mask = _string_missing_mask(df["compound_name"])
        if mask.any() and "molecule_pref_name" in df.columns:
            df.loc[mask, "compound_name"] = df.loc[mask, "molecule_pref_name"].astype(
                "string"
            )
            filled.add("compound_name")

    if "compound_key" in df.columns:
        missing_mask = _string_missing_mask(df["compound_key"])
    else:
        missing_mask = pd.Series(True, index=df.index, dtype="boolean")
        df["compound_key"] = pd.Series(pd.NA, index=df.index, dtype="string")

    if missing_mask.any():
        for candidate in ("parent_molecule_chembl_id", "molecule_chembl_id"):
            if candidate not in df.columns:
                continue
            available_mask = missing_mask & ~_string_missing_mask(df[candidate])
            if available_mask.any():
                df.loc[available_mask, "compound_key"] = df.loc[
                    available_mask, candidate
                ].astype("string")
                filled.add("compound_key")
                missing_mask = _string_missing_mask(df["compound_key"])
            if not missing_mask.any():
                break

    log_value_column_present = "log_value" in df.columns
    log_value_series: pd.Series = (
        pd.to_numeric(df["log_value"], errors="coerce").astype("Float64")
        if log_value_column_present
        else pd.Series(pd.NA, index=df.index, dtype="Float64")
    )
    log_value_filled = False
    pchembl_series: pd.Series | None = None
    if "pchembl_value" in df.columns:
        pchembl_series = pd.to_numeric(df["pchembl_value"], errors="coerce").astype(
            "Float64"
        )
        mask = log_value_series.isna() & pchembl_series.notna()
        if mask.any():
            log_value_series = log_value_series.copy()
            log_value_series.loc[mask] = pchembl_series.loc[mask]
            log_value_filled = True
            filled.add("log_value")

    if log_value_column_present or log_value_series.notna().any():
        # ``log_value`` may originate from CSV exports where the column is typed as
        # ``string``.  Pandas prohibits assigning floats into a ``string`` dtype
        # column which resulted in ``TypeError: Invalid value for dtype 'str'``
        # when we attempted to populate the computed numeric values.
        #
        # By explicitly re-creating the column with the nullable float dtype we
        # guarantee that subsequent partial assignments (for example the
        # pChEMBL-derived fallbacks) work consistently irrespective of the
        # original storage dtype.
        df = df.drop(columns=["log_value"], errors="ignore")
        df["log_value"] = log_value_series.astype("Float64")

    if "standard_value" in df.columns and "standard_units" in df.columns:
        std_value = pd.to_numeric(df["standard_value"], errors="coerce")
        units = df["standard_units"].astype("string").str.strip().str.lower()
        factors = units.map(
            {
                "m": 1.0,
                "mm": 1e-3,
                "um": 1e-6,
                "µm": 1e-6,
                "nm": 1e-9,
                "pm": 1e-12,
            }
        )
        factors = pd.to_numeric(factors, errors="coerce")
        molar = std_value * factors
        computed = pd.Series(pd.NA, index=df.index, dtype="Float64")
        valid = molar.notna() & (molar > 0)
        if valid.any():
            computed_values = -np.log10(molar[valid].astype("float64"))
            computed.loc[valid] = pd.Series(
                computed_values, index=df.index[valid]
            ).astype("Float64")
        mask = log_value_series.isna() & computed.notna()
        if mask.any():
            log_value_series = log_value_series.copy()
            log_value_series.loc[mask] = computed.loc[mask]
            log_value_filled = True
            filled.add("log_value")

    if log_value_column_present or log_value_series.notna().any():
        df = df.drop(columns=["log_value"], errors="ignore")
        df["log_value"] = log_value_series.astype("Float64")

    if (
        not log_value_filled
        and not log_value_column_present
        and not log_value_series.notna().any()
    ):
        filled.discard("log_value")

    bool_defaults = (
        "approx_cited_activity",
        "shuffled_cit",
        "exact_cited_activity",
        "higly_correlated_cit",
        "review_doc",
        "rounded_data_citation",
    )
    for column in bool_defaults:
        if column not in df.columns:
            df[column] = pd.Series(False, index=df.index, dtype="boolean")
            filled.add(column)
            continue

        series = _safe_to_bool(df[column], column)
        try:
            series = series.astype("boolean")
        except (TypeError, ValueError):
            mask = series.isna()
            if mask.any():
                series = series.copy()
                series.loc[mask] = False
                filled.add(column)
            df[column] = series
            continue

        mask = series.isna()
        if mask.any():
            series = series.copy()
            series.loc[mask] = False
            filled.add(column)
        df[column] = series

    string_defaults = {"original_activity_approx", "original_activity_exact"}
    for column in string_defaults:
        if column not in df.columns:
            df[column] = pd.Series(pd.NA, index=df.index, dtype="string")
        else:
            df[column] = df[column].astype("string")

    if "salt_chembl_id" not in df.columns:
        df["salt_chembl_id"] = pd.Series(pd.NA, index=df.index, dtype="string")
        filled.add("salt_chembl_id")
    else:
        mask = _string_missing_mask(df["salt_chembl_id"])
        if mask.any() and "molecule_chembl_id" in df.columns:
            df.loc[mask, "salt_chembl_id"] = df.loc[mask, "molecule_chembl_id"].astype(
                "string"
            )
            filled.add("salt_chembl_id")

    if "nstereo" not in df.columns:
        df["nstereo"] = pd.Series(pd.NA, index=df.index, dtype="Int64")
        filled.add("nstereo")

    return df, filled


def _transform_activity_frame(
    frame: pd.DataFrame,
    *,
    dictionary_root: Path,
    targets_override: Path | None,
) -> pd.DataFrame:
    df, ensured = _ensure_required_input_columns(frame)
    df, identifier_filled = _ensure_compound_key_sources(
        df, dictionary_root=dictionary_root
    )
    df, augmented = _augment_activity_frame(df)
    filled = ensured | identifier_filled | augmented
    missing = _REQUIRED_COLUMNS - set(df.columns)
    if missing:
        available = ", ".join(sorted(df.columns))
        missing_list = ", ".join(sorted(missing))
        raise ActivityExtendedError(
            "Activity export missing required columns: "
            f"{missing_list}. Available columns: {available}"
        )

    df, ensured_filled = _augment_activity_frame(df)
    original_sources, original_identifier_filled = _ensure_compound_key_sources(
        frame, dictionary_root=dictionary_root
    )
    original_filled = _augment_activity_frame(original_sources)[1]
    filled = (
        ensured_filled
        | original_identifier_filled
        | original_filled
        | identifier_filled
    )

    if filled:
        unresolved_columns = sorted(
            column
            for column in filled
            if column not in df.columns
            or (
                df[column].isna().all()
                and column not in _OPTIONAL_EMPTY_BACKFILL_COLUMNS
            )
        )
        resolved_columns = sorted(set(filled) - set(unresolved_columns))

        if resolved_columns:
            logger.info(
                "activity_extended_missing_columns_filled",
                columns=resolved_columns,
            )

        if unresolved_columns:
            logger.warning(
                "activity_extended_missing_columns_unresolved",
                columns=unresolved_columns,
            )

    _resolve_targets_path(dictionary_root, targets_override)

    df = _prepare_unknown_chirality(df)
    df = _apply_multimol_logic(df)
    df = _merge_document_metadata(df, dictionary_root)
    df = _merge_assay_metadata(df, dictionary_root)
    df = _merge_testitem_metadata(df, dictionary_root)
    df = _rename_columns(df)
    df = _compute_citation_flags(df)
    df = _annotate_high_citation(df, dictionary_root)
    df = _extract_activity_properties_flags(df)
    df = _drop_unused_columns(df)
    df = _merge_target_metadata(
        df, dictionary_root=dictionary_root, targets_override=targets_override
    )
    df = _select_and_cast(df)
    return df


def _derive_output_path(input_path: Path) -> Path:
    base_name = helpers.normalise_export_basename(input_path)
    candidate = base_name
    while candidate.startswith("extended."):
        candidate = candidate[len("extended.") :]

    normalised = _normalised_activity_basename(Path(candidate))
    if normalised is not None:
        stamp, _ = normalised
        final_name = f"extended.output.activity_{stamp}.csv"
    else:
        logger.warning("activity_extended_unexpected_basename", basename=base_name)
        normalised_candidate = helpers.normalise_export_basename(Path(candidate))
        final_name = f"extended.{normalised_candidate}"

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
            raise ActivityExtendedError(
                f"Activity export not found: {explicit_input!s}"
            )

    resolved_search_dir = (
        Path(search_dir)
        if search_dir is not None
        else (
            explicit_input.parent
            if explicit_input is not None
            else _current_default_search_dir()
        )
    )
    if not resolved_search_dir.exists():
        raise ActivityExtendedError(
            f"Search directory does not exist: {resolved_search_dir!s}"
        )
    if not resolved_search_dir.is_dir():
        raise ActivityExtendedError(
            f"Search directory is not a directory: {resolved_search_dir!s}"
        )

    input_path = explicit_input or _latest_activity_export(resolved_search_dir)
    dictionary_root = _resolve_dictionary_root(dictionary_dir)

    frame = helpers.read_csv_strict(
        input_path,
        encoding=helpers.ENCODING_FALLBACKS,
        dtype_map=_ACTIVITY_INPUT_SCHEMA,
        na_values=_NA_MARKERS,
    )
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
        non_null = (
            int(processed[column].notna().sum()) if column in processed.columns else 0
        )
        logger.info(
            "activity_extended_non_null_counts",
            column=column,
            non_null=non_null,
            total=len(processed),
        )
    processed = dedupe_final(processed)
    logger.info("activity_extended_saving", path=str(output_path))
    helpers.write_csv(processed, output_path, columns=_FINAL_COLUMN_ORDER)
    logger.info(
        "activity_extended_saved",
        path=str(output_path),
        rows=len(processed),
    )
    return output_path


__all__ = ["process_activity_extended", "ActivityExtendedError"]
