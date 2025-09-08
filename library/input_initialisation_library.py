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
    "exact_cited_activity_samedoc",
    "approx_cited_activity_samedoc",
    "survives_main_steps",
    "num_citations",
    "higly_correlated_cit_samedoc",
    "multimol",
    "publication_date",
    "sort_order.document",
    "document_contains_external_links",
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


def process_activity_table(
    df_activity: pd.DataFrame, dictionary_dir: Path | str
) -> pd.DataFrame:
    """Transform the combined activity table.

    Parameters
    ----------
    df_activity:
        Deduplicated activity dataframe.
    dictionary_dir:
        Directory containing ``targets_type.csv`` and
        ``citation_fraction.csv``.

    Returns
    -------
    pandas.DataFrame
        Processed activity table ready for CSV export. The result mirrors
        the original Power Query logic including citation significance, target
        type flags, IUPHAR class metadata and the ``multifunctional_enzyme``
        indicator.
    """

    logger.info("Processing activity table")
    df = df_activity.copy()

    # --- unknown chirality -------------------------------------------------
    if "nstereo" in df.columns:
        df["unknown_chirality"] = (
            df["nstereo"].astype("Int64").ne(1).fillna(True)
        )
        df.drop(columns=["nstereo"], inplace=True)
    else:
        # If the source table lacks ``nstereo`` values, default to True. This
        # mirrors the Power Query behaviour where missing stereochemistry is
        # treated as unknown.
        df["unknown_chirality"] = pd.Series(True, index=df.index, dtype="boolean")

    # --- multimol assay map ------------------------------------------------
    group_cols = [
        "salt_chembl_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "assay_chembl_id",
        "standard_type",
    ]
    missing = set(group_cols) - set(df.columns)
    if missing:
        raise KeyError(
            "activity table missing columns for multimol grouping: "
            + ", ".join(sorted(missing))
        )

    counts = (
        df.groupby(group_cols, dropna=False)
        .size()
        .rename("Count")
        .reset_index()
    )
    df = df.merge(counts, on=group_cols, how="left")

    mask = (
        df["unknown_chirality"].fillna(True).eq(False)
        & df["multmol_assay"].isna()
        & (df["Count"] > 1)
    )
    assays = set(df.loc[mask, "assay_chembl_id"].dropna().astype(str))
    df["multimol_assay_same"] = df["assay_chembl_id"].isin(assays)

    df["multmol_assay"] = (
        df["multmol_assay"].astype("boolean").fillna(False).astype(bool)
        | df["multimol_assay_same"]
    )
    df.drop(columns=["multimol_assay_same", "Count"], inplace=True)

    # --- rename columns ----------------------------------------------------
    df = df.rename(
        columns={
            "approx_cited_activity": "rounded_data_citation",
            "shuffled_cit": "shuffled_assay",
            "exact_cited_activity": "exact_data_citation",
            "higly_correlated_cit": "higly_correlated_assay",
            "review_doc": "review",
        }
    )

    df = df.rename(
        columns={
            "activity_chembl_id": "activity_id",
            "salt_chembl_id": "saltform_id",
            "molecule_chembl_id": "molecule_id",
            "target_chembl_id": "target_id",
            "assay_chembl_id": "assay_id",
            "document_chembl_id": "document_id",
            "standard_type": "mesurement_type",
            "log_value": "pA_value",
        }
    )

    # --- drop unused columns -----------------------------------------------
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
    df.drop(columns=[c for c in drop_cols if c in df.columns], inplace=True)

    # --- compute citation flags --------------------------------------------
    bool_cols = [
        "exact_data_citation",
        "higly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
    ]
    for col in bool_cols:
        if col in df.columns:
            df[col] = df[col].astype("boolean").fillna(False).astype(bool)

    df["is_citation"] = df[bool_cols].any(axis=1)

    # --- documents with significant citations ------------------------------
    counts_doc = (
        df.groupby("document_id")["is_citation"]
        .agg(n_citation="sum", n_non_citation=lambda s: (~s).sum())
        .reset_index()
    )
    counts_doc["N"] = counts_doc["n_citation"] + counts_doc["n_non_citation"]
    counts_doc = counts_doc[
        (counts_doc["n_citation"] > 0) & (counts_doc["n_non_citation"] > 0)
    ]

    cf_path = Path(dictionary_dir) / "citation_fraction.csv"
    cf = pd.read_csv(cf_path)

    counts_doc = counts_doc.merge(cf[["N", "K_min_significant"]], on="N", how="left")
    counts_doc["high_citation_rate"] = counts_doc["K_min_significant"].notna() & (
        counts_doc["n_citation"] >= counts_doc["K_min_significant"]
    )

    df = df.merge(
        counts_doc[["document_id", "high_citation_rate"]],
        on="document_id",
        how="left",
    )
    df["high_citation_rate"] = df["high_citation_rate"].fillna(False)

    # --- target types ------------------------------------------------------
    targets_path = Path(dictionary_dir) / "targets_type.csv"
    targets = pd.read_csv(
        targets_path,
        dtype={
            "chembl_id": "string",
            "IUPHAR_class": "string",
            "IUPHAR_subclass": "string",
            "type": "string",
        },
    )

    df = df.merge(
        targets[["chembl_id", "IUPHAR_class", "IUPHAR_subclass", "type"]],
        how="left",
        left_on="target_id",
        right_on="chembl_id",
    )
    mapping = {
        "Multicellular organism": False,
        "Viruses": True,
        "Unicellular organism": True,
    }
    df["unicellular_organism"] = (
        df["type"].map(mapping).fillna(False).astype(bool)
    )
    df["multifunctional_enzyme"] = df["IUPHAR_subclass"].eq("Multifunctional")
    df.drop(columns=["chembl_id", "type"], inplace=True)

    # --- final ordering ----------------------------------------------------
    final_cols = [
        "activity_id",
        "saltform_id",
        "molecule_id",
        "target_id",
        "assay_id",
        "document_id",
        "bao_endpoint",
        "mesurement_type",
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
    ]

    missing_final = set(final_cols) - set(df.columns)
    if missing_final:
        raise KeyError(
            "activity table missing expected columns: "
            + ", ".join(sorted(missing_final))
        )
    df = df[final_cols]

    # Ensure boolean dtype where appropriate
    bool_cols_final = [
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
    for col in bool_cols_final:
        df[col] = df[col].astype("boolean")

    return df


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
    dictionary_dir: Path | str | None = None,
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
            if dictionary_dir is not None:
                df = process_activity_table(df, dictionary_dir)
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
