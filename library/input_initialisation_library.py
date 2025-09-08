"""Utilities for combining ChEMBL initialisation input tables.

This module provides helpers to read source Excel workbooks, softly cast
common column types, merge entity tables and persist the final CSV files.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Literal

import pandas as pd

logger = logging.getLogger(__name__)

EntityName = Literal["activity", "assay", "document", "target", "testitem"]

# Mapping of sheet names to entity identifiers
SAME_DOC_SHEETS: dict[str, EntityName] = {
    "assay_step5_same_doc": "assay",
    "molecule_step5_same_doc": "testitem",
    "target_step5_same_doc": "target",
    "document_step5_same_doc": "document",
    "activities_step5_same_doc": "activity",
}

ALL_DOC_SHEETS: dict[str, EntityName] = {
    "assay_step5": "assay",
    "molecule_step5": "testitem",
    "target_step5": "target",
    "document_step5": "document",
    "activities_step5": "activity",
}

# Columns that should be removed from the combined activity table
ACTIVITY_DROP_COLS: set[str] = {
    "n_stereocenters",
    "unspecified_chirality_mol",
    "shuffled_target_assay",
    "approx_cited_activity",
    "exact_cited_activity_samedoc",
    "approx_cited_activity_samedoc",
    "survives_main_steps",
    "num_citations",
    "higly_correlated_cit_samedoc",
    "higly_correlated_cit",
    "shuffled_cit",
    "exact_cited_activity",
    "multmol_assay",
    "multimol",
    "publication_date",
    "sort_order.document",
    "document_contains_external_links",
    "review_doc",
    "year",
    "month",
    "day",
    "Citation",
    "Column1",
    "Column2",
    "Column3",
    "Column4",
    "Column5",
    "Column6",
    "Column7",
    "Column8",
    "Column9",
    "Column10",
    "Column11",
    "Column12",
    "Column13",
    "Column14",
    "Column15",
    "Column16",
    "Column17",
    "Column18",
    "Column19",
    "Column20",
    "Column21",
    "Column22",
    "Column23",
    "Column24",
    "Column25",
    "Column26",
    "Column27",
    "Column28",
    "Column29",
    "Column30",
    "Column31",
    "Column32",
    "Column33",
    "Column34",
    "Column35",
    "Column36",
    "Column37",
    "Column38",
    "Column39",
    "Column40",
    "Column41",
    "Column42",
    "Column43",
    "Column44",
}

# Type hints for soft conversions
STRING_COLS: set[str] = {
    "assay_chembl_id",
    "document_chembl_id",
    "target_chembl_id",
    "salt_chembl_id",
    "molecule_chembl_id",
}

INT_COLS: set[str] = {
    "isoform",
    "version",
    "acts_per_assay_step5",
    "volume",
    "last_page",
    "pubmed_id",
    "classification",
    "year",
    "month",
    "day",
    "citation",
    "nActivity",
    "nstereo",
    "n_stereocenters",
}

BOOL_COLS: set[str] = {
    "shuffled_target_assay",
    "shuffled_cit",
    "higly_correlated_cit",
    "cited_assay_corr",
    "error_assay_corr",
    "document_contains_external_links",
    "is_experimental_doc",
}

FLOAT_COLS: set[str] = {
    "mw_freebase",
    "standard_value",
    "log_value",
    "citation_fraction",
}

DATE_COLS: set[str] = {
    "publication_date",
}


def _ensure_openpyxl() -> None:
    """Ensure a supported ``openpyxl`` version is available.

    Raises
    ------
    RuntimeError
        If ``openpyxl`` is missing or too old.
    """
    try:
        import openpyxl  # type: ignore
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError("openpyxl>=3.1 is required to read Excel files") from exc
    version = tuple(int(v) for v in openpyxl.__version__.split(".")[:3])
    if version < (3, 1, 0):  # pragma: no cover - defensive
        raise RuntimeError(f"openpyxl>=3.1 is required, found {openpyxl.__version__}")


def _read_sheet(
    path: Path, sheet: str, *, header_promotion: bool = False
) -> pd.DataFrame:
    """Read ``sheet`` from an Excel file.

    Parameters
    ----------
    path:
        Excel workbook location.
    sheet:
        Sheet name to load.
    header_promotion:
        If ``True`` the first row is treated as header regardless of Excel's
        stored header information.

    Raises
    ------
    ValueError
        If the sheet is missing.
    """
    _ensure_openpyxl()
    try:
        if header_promotion:
            raw = pd.read_excel(path, sheet_name=sheet, header=None, engine="openpyxl")
            if raw.empty:
                return pd.DataFrame()
            header = raw.iloc[0].astype(str).tolist()
            df = raw.iloc[1:].reset_index(drop=True)
            df.columns = header
            return df
        return pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    except ValueError as exc:  # missing sheet
        raise ValueError(f"sheet '{sheet}' not found in {path}") from exc


def load_same_doc(xlsx_path: Path) -> Dict[EntityName, pd.DataFrame]:
    """Load entity tables from ``ChEMBL_same_document_20_05.xlsx``.

    Parameters
    ----------
    xlsx_path:
        Path to the Excel workbook.

    Returns
    -------
    dict
        Mapping of entity names to DataFrames.
    """
    tables: Dict[EntityName, pd.DataFrame] = {}
    for sheet, entity in SAME_DOC_SHEETS.items():
        logger.info("Reading sheet '%s' from %s", sheet, xlsx_path)
        tables[entity] = _read_sheet(xlsx_path, sheet)
    return tables


def load_all_doc(xlsx_path: Path) -> Dict[EntityName, pd.DataFrame]:
    """Load entity tables from ``ChEMBL_all_10_05_step5.xlsx``.

    The ``activities_step5`` sheet requires promotion of the first row to
    header.
    """
    tables: Dict[EntityName, pd.DataFrame] = {}
    for sheet, entity in ALL_DOC_SHEETS.items():
        logger.info("Reading sheet '%s' from %s", sheet, xlsx_path)
        promote = sheet == "activities_step5"
        tables[entity] = _read_sheet(xlsx_path, sheet, header_promotion=promote)
    return tables


def _safe_to_int(series: pd.Series, col: str) -> pd.Series:
    try:
        return pd.to_numeric(series, errors="raise").astype("Int64")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("Failed to convert column '%s' to Int64: %s", col, exc)
        return series.astype("string")


def _safe_to_float(series: pd.Series, col: str) -> pd.Series:
    try:
        return pd.to_numeric(series, errors="raise").astype("float64")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("Failed to convert column '%s' to float: %s", col, exc)
        return series.astype("string")


def _safe_to_datetime(series: pd.Series, col: str) -> pd.Series:
    try:
        return pd.to_datetime(series, errors="raise")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("Failed to convert column '%s' to datetime: %s", col, exc)
        return series.astype("string")


def _safe_to_bool(series: pd.Series, col: str) -> pd.Series:
    def mapper(value: Any) -> object:
        if pd.isna(value):
            return pd.NA
        if isinstance(value, str):
            value = value.strip().lower()
        if value in {True, 1, "1", "true", "t"}:
            return True
        if value in {False, 0, "0", "false", "f"}:
            return False
        raise ValueError(f"invalid boolean value: {value}")

    try:
        mapped = series.map(mapper)
        return mapped.astype("boolean")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("Failed to convert column '%s' to boolean: %s", col, exc)
        return series.astype("string")


def unify_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Softly cast frequent columns to expected dtypes.

    Unrecognised values fall back to ``string`` dtype while emitting a warning.
    """
    result = df.copy()
    for col in STRING_COLS & set(result.columns):
        result[col] = result[col].astype("string")
    for col in INT_COLS & set(result.columns):
        result[col] = _safe_to_int(result[col], col)
    for col in BOOL_COLS & set(result.columns):
        result[col] = _safe_to_bool(result[col], col)
    for col in FLOAT_COLS & set(result.columns):
        result[col] = _safe_to_float(result[col], col)
    for col in DATE_COLS & set(result.columns):
        result[col] = _safe_to_datetime(result[col], col)
    return result


def append_entities(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Vertically concatenate ``df_a`` and ``df_b`` and remove duplicates."""
    combined = pd.concat([df_a, df_b], ignore_index=True, sort=False)
    return combined.drop_duplicates(keep="first")


def build_combined_tables(
    same: Dict[EntityName, pd.DataFrame],
    all_: Dict[EntityName, pd.DataFrame],
) -> Dict[EntityName, pd.DataFrame]:
    """Combine entity tables from ``same`` and ``all_`` sources."""
    combined: Dict[EntityName, pd.DataFrame] = {}
    for entity in SAME_DOC_SHEETS.values():
        df_same = unify_dtypes(same[entity])
        df_all = unify_dtypes(all_[entity])
        rows_same = len(df_same)
        rows_all = len(df_all)
        if entity == "activity":
            concat = pd.concat([df_same, df_all], ignore_index=True, sort=False)
            concat = concat.drop(
                columns=list(ACTIVITY_DROP_COLS & set(concat.columns)),
                errors="ignore",
            )
            df = concat.drop_duplicates(keep="first")
            rows_concat = len(concat)
            rows_after = len(df)
        else:
            df = append_entities(df_same, df_all)
            rows_concat = rows_same + rows_all
            rows_after = len(df)
        logger.info(
            "Entity %s: rows_same=%d rows_all=%d rows_concat=%d rows_after_dedup=%d",
            entity,
            rows_same,
            rows_all,
            rows_concat,
            rows_after,
        )
        if entity == "activity":
            remaining = ACTIVITY_DROP_COLS & set(df.columns)
            if remaining:
                raise ValueError(
                    f"activity table still contains columns: {', '.join(sorted(remaining))}"
                )
        combined[entity] = df
    return combined


def save_tables(
    tables: Dict[EntityName, pd.DataFrame],
    out_dir: Path,
    fmt: str = "csv",
) -> Dict[EntityName, Path]:
    """Persist combined tables to ``out_dir``.

    Parameters
    ----------
    tables:
        Mapping of entity names to DataFrames.
    out_dir:
        Destination directory. Created if missing.
    fmt:
        Output format. Only ``"csv"`` is supported.

    Returns
    -------
    dict
        Mapping of entity names to written file paths.
    """
    if fmt != "csv":
        raise ValueError("only csv output is supported")
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[EntityName, Path] = {}
    for entity, df in tables.items():
        path = out_dir / f"{entity}.csv"
        df.to_csv(path, index=False, encoding="utf-8", na_rep="")
        logger.info("Wrote %d rows to %s", len(df), path)
        paths[entity] = path
    return paths
