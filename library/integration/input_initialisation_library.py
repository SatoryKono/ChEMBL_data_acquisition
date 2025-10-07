"""Utilities for combining ChEMBL initialisation input tables.

This module provides helpers to read source Excel workbooks, softly cast
common column types, merge entity tables and persist the final CSV files.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, Literal

try:
    from dataclasses import dataclass
except ImportError as exc:  # pragma: no cover - Python <3.7
    raise ImportError(
        "dataclasses module is required; upgrade to Python 3.7 or newer"
    ) from exc

import pandas as pd

from ..common.log import logger
from ..config import Config
from ..io.writers import write_csv
from ..pipelines.target import organism_classification

EntityName = Literal[
    "activity",
    "assay",
    "document",
    "target",
    "testitem",
    "pairs",
    "pairs_same_document",
]

# Mapping of entity names to their corresponding dataframes.  The dictionary
# may contain additional keys produced during processing, hence the generic
# ``str`` key type.
TableDict = dict[str, pd.DataFrame]


# Mapping of sheet names to entity identifiers
SAME_DOC_SHEETS: dict[str, EntityName] = {
    "assay_step5_same_doc": "assay",
    "molecule_step5_same_doc": "testitem",
    "target_step5_same_doc": "target",
    "document_step5_same_doc": "document",
    "activities_step5_same": "activity",
    "pairs_same_doc": "pairs_same_document",
}

ALL_DOC_SHEETS: dict[str, EntityName] = {
    "assay_step5": "assay",
    "molecule_step5": "testitem",
    "target_step5": "target",
    "document_step5": "document",
    "activities_step5": "activity",
    "step5_pairs": "pairs",
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

TARGET_TYPE_USECOLS: frozenset[str] = frozenset(
    {
        "target_chembl_id",
        "target_sort_order",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "gene_index",
        "taxon_index",
        "multifunctional_enzyme",
        "genus",
        "superkingdom",
        "phylum",
    }
)

TARGET_TYPE_OPTIONAL_COLS: frozenset[str] = frozenset({"lineage_class"})

TARGET_TYPE_DTYPES: dict[str, str] = {
    "target_chembl_id": "string",
    "target_sort_order": "string",
    "IUPHAR_class": "string",
    "IUPHAR_subclass": "string",
    "gene_index": "string",
    "taxon_index": "string",
    "multifunctional_enzyme": "string",
    "genus": "string",
    "superkingdom": "string",
    "phylum": "string",
    "lineage_class": "string",
}


def get_percentage(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """Calculate percentage distribution for a table with a ``Filtered`` column.

    Parameters
    ----------
    df:
        DataFrame containing a ``Filtered`` column.
    table_name:
        Logical name of the table for error reporting.

    Returns
    -------
    pandas.DataFrame
        Table with columns ``Filtered``, ``Count`` and ``Percentage, %``
        including a summary ``Total`` row.

    Raises
    ------
    KeyError
        If the ``Filtered`` column is absent.
    """

    if "Filtered" not in df.columns:
        available = ", ".join(df.columns)
        raise KeyError(
            f"table '{table_name}' missing column 'Filtered'; available: {available}"
        )

    counts = df.groupby("Filtered", dropna=False).size().rename("Count").reset_index()
    total = int(counts["Count"].sum())

    if total > 0:
        counts["fraction"] = counts["Count"] / total * 100
    else:
        counts["fraction"] = 0.0

    total_row = pd.DataFrame(
        {"Filtered": ["Total"], "Count": [total], "fraction": [100.0 if total else 0.0]}
    )
    counts = pd.concat([counts, total_row], ignore_index=True)

    def _round(val: float) -> float:
        if val == 100 or val == total:
            return float(total if total else 0)
        if val > 10:
            return round(val, 1)
        if val > 1:
            return round(val, 2)
        if val > 0.1:
            return round(val, 3)
        return round(val, 4)

    counts["Percentage, %"] = counts["fraction"].apply(_round)
    return counts.drop(columns=["fraction"])


def add_percentage(
    statistics: pd.DataFrame, percent_df: pd.DataFrame, table_name: str
) -> pd.DataFrame:
    """Merge percentage information into aggregated statistics.

    Parameters
    ----------
    statistics:
        Aggregated metrics per ``Filtered`` value including a ``Total`` row.
    percent_df:
        Output of :func:`get_percentage` for the same table.
    table_name:
        Entity name used to prefix the percentage column.

    Returns
    -------
    pandas.DataFrame
        ``percent_df`` enriched with metric columns and renamed percentage
        column ``{table_name}.Percentage, %``.
    """

    merged = percent_df.merge(statistics, on="Filtered", how="left")
    col_name = f"{table_name}.Percentage, %"
    merged.rename(columns={"Percentage, %": col_name}, inplace=True)

    metric_cols = [c for c in statistics.columns if c != "Filtered"]
    ordered = ["Filtered", "Count", *metric_cols, col_name]
    return merged[ordered]


def compute_status_statistics(df: pd.DataFrame, table_name: str) -> pd.DataFrame:
    """Prepare statistics with percentage distribution.

    Parameters
    ----------
    df:
        DataFrame containing ``Filtered.new`` and metric columns.
    table_name:
        Entity name used for percentage column prefix.

    Returns
    -------
    pandas.DataFrame
        Aggregated table grouped by ``Filtered`` with counts and percentage
        information appended.
    """

    if "Filtered.new" not in df.columns:
        available = ", ".join(df.columns)
        raise KeyError(
            f"table '{table_name}' missing column 'Filtered.new'; "
            f"available: {available}"
        )

    df_tmp = df.rename(columns={"Filtered.new": "Filtered"}).copy()
    metric_cols = [
        c
        for c in df_tmp.columns
        if c != "Filtered" and pd.api.types.is_numeric_dtype(df_tmp[c])
    ]

    grouped = df_tmp.groupby("Filtered", dropna=False)[metric_cols].sum().reset_index()
    totals = {col: grouped[col].sum() for col in metric_cols}
    grouped = pd.concat(
        [grouped, pd.DataFrame([{"Filtered": "Total", **totals}])],
        ignore_index=True,
    )

    percent_df = get_percentage(df_tmp, table_name)
    return add_percentage(grouped, percent_df, table_name)


def process_activity_table(
    df_activity: pd.DataFrame,
    dictionary_dir: Path | str,
    targets_csv: Path | str | None = None,
) -> pd.DataFrame:
    """Transform the combined activity table.

    Parameters
    ----------
    df_activity:
        Deduplicated activity dataframe.
    dictionary_dir:

        Directory containing ``_Curation/citation_fraction.csv`` and
        ``targets_type.csv`` in a ``_target`` subdirectory.

    targets_csv:
        Optional explicit path to ``targets_type.csv``. When provided, the file
        is loaded from this location instead of searching within
        ``dictionary_dir``.

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
        df["unknown_chirality"] = df["nstereo"].astype("Int64").ne(1).fillna(True)
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

    counts = df.groupby(group_cols, dropna=False).size().rename("Count").reset_index()
    df = df.merge(counts, on=group_cols, how="left")

    mask = (
        df["unknown_chirality"].fillna(True).eq(False)
        & df["multmol_assay"].isna()
        & (df["Count"] > 1)
    )
    assays = set(df.loc[mask, "assay_chembl_id"].dropna().astype(str))
    df["multimol_assay_same"] = df["assay_chembl_id"].isin(assays)

    multmol_assay_series = _safe_to_bool(df["multmol_assay"], "multmol_assay").fillna(
        False
    )
    df["multmol_assay"] = _safe_to_bool(
        multmol_assay_series | df["multimol_assay_same"], "multmol_assay"
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
            "activity_chembl_id": "activity_chembl_id",
            "salt_chembl_id": "saltform_id",
            # keep CHEMBL identifiers for related entities
            "standard_type": "standard_type",
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
            df[col] = _safe_to_bool(df[col], col).fillna(False)

    df["is_citation"] = df[bool_cols].any(axis=1)

    # --- documents with significant citations ------------------------------
    counts_doc = (
        df.groupby("document_chembl_id")["is_citation"]
        .agg(n_citation="sum", n_non_citation=lambda s: (~s).sum())
        .reset_index()
    )
    counts_doc["N"] = counts_doc["n_citation"] + counts_doc["n_non_citation"]
    counts_doc = counts_doc[
        (counts_doc["n_citation"] > 0) & (counts_doc["n_non_citation"] > 0)
    ]

    cf_path = Path(dictionary_dir) / "_Curation" / "citation_fraction.csv"
    cf = pd.read_csv(cf_path)

    counts_doc = counts_doc.merge(cf[["N", "K_min_significant"]], on="N", how="left")
    counts_doc["high_citation_rate"] = counts_doc["K_min_significant"].notna() & (
        counts_doc["n_citation"] >= counts_doc["K_min_significant"]
    )

    df = df.merge(
        counts_doc[["document_chembl_id", "high_citation_rate"]],
        on="document_chembl_id",
        how="left",
    )
    df["high_citation_rate"] = _safe_to_bool(
        df["high_citation_rate"], "high_citation_rate"
    ).fillna(False)

    # --- target types ------------------------------------------------------
    if targets_csv is not None:
        targets_path = Path(targets_csv)
    else:
        targets_path = Path(dictionary_dir) / "targets_type.csv"
        if not targets_path.exists():
            targets_path = Path(dictionary_dir) / "_target" / "targets_type.csv"
        if not targets_path.exists():
            legacy_path = Path(dictionary_dir) / "_Target" / "targets_type.csv"
            if legacy_path.exists():
                targets_path = legacy_path
            else:
                expected_locations = [
                    Path(dictionary_dir) / "targets_type.csv",
                    Path(dictionary_dir) / "_target" / "targets_type.csv",
                ]
                msg = (
                    "targets_type.csv not found in the provided dictionary directory. "
                    "Expected at either "
                    + " or ".join(f"'{path}'" for path in expected_locations)
                    + "."
                )
                raise FileNotFoundError(msg)

    targets = pd.read_csv(
        targets_path,
        usecols=[
            "target_chembl_id",
            "IUPHAR_class",
            "IUPHAR_subclass",
            "gene_index",
            "taxon_index",
            "target_sort_order",
            "multifunctional_enzyme",
            "genus",
            "superkingdom",
            "phylum",
            "taxon_id",
        ],
        dtype={
            "target_chembl_id": "string",
            "IUPHAR_class": "string",
            "IUPHAR_subclass": "string",
            "taxon_index": "string",
            "gene_index": "string",
            "target_sort_order": "string",
            "multifunctional_enzyme": "string",
            "genus": "string",
            "superkingdom": "string",
            "phylum": "string",
            "taxon_id": "string",
        },
    )

    df = df.merge(
        targets[
            [
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
            ]
        ],
        how="left",
        on="target_chembl_id",
    )
    # Drop any duplicated columns to avoid unexpected suffixed names when the
    # input dataframe already contains some of these fields.
    df = df.loc[:, ~df.columns.duplicated()]

    if "organism_cellularity" not in df.columns:
        df["organism_cellularity"] = pd.Series(dtype="string")

    df["multifunctional_enzyme"] = _safe_to_bool(
        df["multifunctional_enzyme"], "multifunctional_enzyme"
    )

    df = organism_classification.add_cellularity_smart(
        df,
        genus_col="genus",
        superkingdom_col="superkingdom",
        phylum_col="phylum",
        output_col="organism_cellularity",
    )

    unicellular_labels = {
        organism_classification.TYPE_UNICELLULAR,
        organism_classification.TYPE_VIRAL,
    }
    df["unicellular_organism"] = (
        df["organism_cellularity"].astype("string").isin(unicellular_labels)
    ).astype("boolean")

    df.drop(
        columns=[
            "genus",
            "superkingdom",
            "phylum",
            "taxon_id",
            "organism_cellularity",
        ],
        inplace=True,
    )

    # --- final ordering ----------------------------------------------------
    final_cols = [
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
    ]

    # Remove any accidental duplicates while preserving the declared order. This
    # guards against repeated column names that would otherwise lead to
    # duplicated output columns downstream.
    final_cols = list(dict.fromkeys(final_cols))

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
        df[col] = _safe_to_bool(df[col], col)

    return df


def _ensure_openpyxl() -> None:
    """Ensure a supported ``openpyxl`` version is available.

    Raises
    ------
    RuntimeError
        If ``openpyxl`` is missing or too old.

    """
    try:
        import openpyxl
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError("openpyxl>=3.1 is required to read Excel files") from exc
    version = tuple(int(v) for v in openpyxl.__version__.split(".")[:3])
    if version < (3, 1, 0):  # pragma: no cover - defensive
        raise RuntimeError(f"openpyxl>=3.1 is required, found {openpyxl.__version__}")


def _read_sheet(
    path: Path, sheet: str, *, header_promotion: bool = False
) -> pd.DataFrame:
    """Read ``sheet`` from an Excel file.

    Duplicate column names are dropped with a warning.

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
            df.columns = pd.Index(header)
        else:
            df = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")

        if df.columns.has_duplicates:
            logger.warning(
                "duplicate columns dropped",
                sheet=sheet,
                duplicates=df.columns[df.columns.duplicated()].tolist(),
            )
            df = df.loc[:, ~df.columns.duplicated()]
        return df
    except ValueError as exc:  # missing sheet
        raise ValueError(f"sheet '{sheet}' not found in {path}") from exc


def load_same_doc(xlsx_path: Path) -> TableDict:
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
    tables: TableDict = {}
    for sheet, entity in SAME_DOC_SHEETS.items():
        logger.info("sheet_read", sheet=sheet, path=str(xlsx_path))
        tables[entity] = _read_sheet(xlsx_path, sheet)
    return tables


def load_all_doc(xlsx_path: Path) -> TableDict:
    """Load entity tables from ``ChEMBL_all_10_05_step5.xlsx``.

    The ``activities_step5`` sheet requires promotion of the first row to
    header.
    """
    tables: TableDict = {}
    for sheet, entity in ALL_DOC_SHEETS.items():
        logger.info("sheet_read", sheet=sheet, path=str(xlsx_path))
        promote = sheet == "activities_step5"
        tables[entity] = _read_sheet(xlsx_path, sheet, header_promotion=promote)
    return tables


def _safe_to_int(series: pd.Series, col: str) -> pd.Series:
    """Safely cast ``series`` to ``Int64``.

    Parameters
    ----------
    series:
        Column values to convert.
    col:
        Column name used in warning messages.

    Returns
    -------
    pandas.Series
        ``Int64`` typed series or ``string`` on failure.
    """

    if not isinstance(series, pd.Series):
        raise TypeError(f"column '{col}' has duplicate entries; expected a Series")

    try:
        return pd.to_numeric(series, errors="raise").astype("Int64")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("dtype_int_conversion_failed", column=col, error=str(exc))
        return series.astype("string")


def _safe_to_float(series: pd.Series, col: str) -> pd.Series:
    """Safely cast ``series`` to ``float64``.

    Parameters
    ----------
    series:
        Column values to convert.
    col:
        Column name used in warning messages.

    Returns
    -------
    pandas.Series
        ``float64`` typed series or ``string`` on failure.
    """

    if not isinstance(series, pd.Series):
        raise TypeError(f"column '{col}' has duplicate entries; expected a Series")

    try:
        return pd.to_numeric(series, errors="raise").astype("float64")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("dtype_float_conversion_failed", column=col, error=str(exc))
        return series.astype("string")


def _safe_to_datetime(series: pd.Series, col: str) -> pd.Series:
    """Safely cast ``series`` to pandas ``datetime64``.

    Parameters
    ----------
    series:
        Column values to convert.
    col:
        Column name used in warning messages.

    Returns
    -------
    pandas.Series
        ``datetime64`` typed series or ``string`` on failure.
    """

    if not isinstance(series, pd.Series):
        raise TypeError(f"column '{col}' has duplicate entries; expected a Series")

    try:
        return pd.to_datetime(series, errors="raise")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("dtype_datetime_conversion_failed", column=col, error=str(exc))
        return series.astype("string")


def _safe_to_bool(series: pd.Series, col: str) -> pd.Series:
    """Safely cast ``series`` to pandas ``boolean``.

    Parameters
    ----------
    series:
        Column values to convert.
    col:
        Column name used in warning messages.

    Returns
    -------
    pandas.Series
        ``boolean`` typed series or ``string`` on failure.
    """

    if not isinstance(series, pd.Series):
        raise TypeError(f"column '{col}' has duplicate entries; expected a Series")

    def mapper(value: Any) -> object:
        if pd.isna(value):
            return pd.NA
        if isinstance(value, str):
            value = value.strip().lower()

        if value in (True, 1, "1", "true", "t"):
            return True
        if value in (False, 0, "0", "false", "f"):
            return False
        raise ValueError(f"invalid boolean value: {value}")

    try:
        mapped = series.map(mapper)
        return mapped.astype("boolean")
    except Exception as exc:  # pragma: no cover - rare
        logger.warning("dtype_bool_conversion_failed", column=col, error=str(exc))
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


def generate_pair_entity_tables(
    tables: TableDict, pair_keys: Mapping[str, str]
) -> TableDict:
    """Generate entity subtables for each pair type.

    For every pair table listed in ``pair_keys`` the function extracts unique
    activity identifiers, filters the ``activity`` table accordingly and
    produces corresponding ``assay``, ``document``, ``target`` and ``testitem``
    subtables. New tables are appended to a copy of ``tables`` using the
    provided suffixes.

    Args:
        tables: Mapping of entity names to their dataframes. Must include the
            ``activity`` table and the pair tables referenced in ``pair_keys``.
        pair_keys: Mapping of pair table keys to suffixes used for naming the
            resulting subtables. For example ``{"pairs_independent": "ind"}``
            would create tables named ``activity_ind``, ``assay_ind`` etc.

    Returns:
        TableDict: Copy of ``tables`` extended with newly generated entity
        tables.

    Warns:
        logger.warning: If required tables or columns are missing. In such
        cases the corresponding subtables are skipped.
    """

    result: TableDict = {**tables}

    activity_df = tables.get("activity")
    if activity_df is None:
        logger.warning("pair_activity_table_missing")
        return result
    if "activity_chembl_id" not in activity_df.columns:
        logger.warning("pair_activity_id_column_missing")
        return result

    # Columns linking activities to related entities use CHEMBL identifiers.
    # This mapping specifies the identifier column for each entity so that pair
    # table generation can correctly retrieve the associated rows.
    # Mapping of entity names to their identifier columns. Test items now use
    # the canonical ``molecule_chembl_id`` to ensure consistent naming across
    # the codebase.
    entity_cols: dict[str, str] = {
        "assay": "assay_chembl_id",
        "document": "document_chembl_id",
        "target": "target_chembl_id",
        "testitem": "molecule_chembl_id",
    }

    for pair_key, suffix in pair_keys.items():
        pairs_df = tables.get(pair_key)
        if pairs_df is None:
            logger.warning("pair_table_missing", pair_table=pair_key)
            continue

        # Harmonise activity identifier column names to ensure downstream
        # processing can rely on a consistent schema regardless of the source
        # workbook's naming conventions.
        pairs_df = normalize_pair_columns(pairs_df)
        result[pair_key] = pairs_df

        required = {"activity_chembl_id1", "activity_chembl_id2"}
        missing = required - set(pairs_df.columns)
        if missing:
            logger.warning(
                "pair_table_missing_columns",
                pair_table=pair_key,
                columns=sorted(missing),
            )
            continue

        activity_ids = pd.unique(
            pd.concat(
                [pairs_df["activity_chembl_id1"], pairs_df["activity_chembl_id2"]],
                ignore_index=True,
            ).dropna()
        )

        filtered_activity = activity_df[
            activity_df["activity_chembl_id"].isin(activity_ids)
        ].copy()
        result[f"activity_{suffix}"] = filtered_activity

        for entity, id_col in entity_cols.items():
            if id_col not in filtered_activity.columns:
                logger.warning(
                    "pair_activity_column_missing",
                    column=id_col,
                    entity=entity,
                    suffix=suffix,
                )
                continue

            ids = pd.unique(filtered_activity[id_col].dropna())
            entity_df = tables.get(entity)
            if entity_df is None:
                logger.warning("pair_entity_table_missing", entity=entity)
                continue
            if id_col not in entity_df.columns:
                logger.warning(
                    "pair_entity_column_missing", entity=entity, column=id_col
                )
                continue
            result[f"{entity}_{suffix}"] = entity_df[entity_df[id_col].isin(ids)].copy()

    return result


def add_pair_metric_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add independent/non-independent metric flags to pair tables.

    Parameters
    ----------
    df:
        DataFrame containing pair information. Expected to include the
        ``INDEPENDENT`` and ``standard_type`` columns.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with four additional integer columns indicating the
        count of independent or non-independent measurements for ``IC50`` and
        ``Ki`` values.

    """
    result = df.copy()
    indep_col, type_col = "INDEPENDENT", "standard_type"
    if {indep_col, type_col}.issubset(result.columns):
        # Convert flag column to pandas boolean and fill missing with ``False``
        independent = _safe_to_bool(result[indep_col], indep_col).fillna(False)
        mtype = result[type_col].astype("string")
        # Derive measurement-specific indicators as integer columns
        result["non_independent_IC50"] = (~independent & (mtype == "IC50")).astype(
            "int64"
        )
        result["independent_IC50"] = (independent & (mtype == "IC50")).astype("int64")
        result["non_independent_Ki"] = (~independent & (mtype == "Ki")).astype("int64")
        result["independent_Ki"] = (independent & (mtype == "Ki")).astype("int64")
    else:
        # Ensure metric columns exist even if source fields are missing
        for col in (
            "non_independent_IC50",
            "independent_IC50",
            "non_independent_Ki",
            "independent_Ki",
        ):
            result[col] = 0
        missing = [c for c in (indep_col, type_col) if c not in result.columns]
        logger.warning("pair_metric_columns_missing", columns=missing)
    return result


def build_combined_tables(
    same: TableDict,
    all_: TableDict,
    dictionary_dir: Path | str | None = None,
    *,
    targets_type_csv: Path | str | None = None,
) -> TableDict:
    """Combine entity tables from ``same`` and ``all_`` sources.

    Deduplication uses specific ChEMBL identifier columns for each entity
    rather than relying on column order.

    Parameters
    ----------
    same:
        Tables generated from the "same document" workbook.
    all_:
        Tables generated from the "all document" workbook.
    dictionary_dir:
        Optional directory containing auxiliary CSV dictionaries used during
        processing.
    targets_type_csv:
        Optional path to ``targets_type.csv``. If relative, it is resolved
        against ``dictionary_dir``.

    Returns
    -------
    TableDict
        Mapping of entity names to combined tables.

    """
    combined: TableDict = {}
    dict_dir = Path(dictionary_dir) if dictionary_dir is not None else None
    targets_path = Path(targets_type_csv) if targets_type_csv is not None else None

    if (
        targets_path is not None
        and not targets_path.is_absolute()
        and dict_dir is not None
    ):
        targets_path = dict_dir / targets_path.name

    ID_COLS: dict[str, str | list[str]] = {
        "assay": "assay_chembl_id",
        "document": "document_chembl_id",
        "target": "target_chembl_id",
        "testitem": ["salt_chembl_id", "saltform_id"],
        "activity": "activity_chembl_id",
    }

    # --- regular entities -------------------------------------------------
    regular_entities: tuple[EntityName, ...] = (
        "assay",
        "document",
        "target",
        "testitem",
    )
    for entity in regular_entities:
        df_same = unify_dtypes(same[entity])
        df_all = unify_dtypes(all_[entity])
        df_tmp = append_entities(df_same, df_all)
        if df_tmp.shape[1]:
            id_col: str | list[str] | None = ID_COLS[entity]
            if isinstance(id_col, list):
                id_col = next((c for c in id_col if c in df_tmp.columns), None)
            df = (
                df_tmp.drop_duplicates(subset=id_col, keep="first")
                if id_col
                else df_tmp
            )
        else:
            df = df_tmp
        logger.info(
            "entity_rows_summary",
            entity=entity,
            rows_same=len(df_same),
            rows_all=len(df_all),
            rows_after_dedup=len(df),
        )
        combined[entity] = df

    # --- activity --------------------------------------------------------
    df_same_act = unify_dtypes(same["activity"])
    df_same_act = df_same_act.loc[:, ~df_same_act.columns.duplicated()]
    df_all_act = unify_dtypes(all_["activity"])

    # Drop duplicate columns from each DataFrame before concatenation. Pandas
    # ``concat`` requires unique column names across inputs.
    same_dups = df_same_act.columns[df_same_act.columns.duplicated()].tolist()
    if same_dups:
        logger.info(
            "activity_duplicate_columns_removed", source="same", columns=same_dups
        )
        df_same_act = df_same_act.loc[:, ~df_same_act.columns.duplicated()]
    all_dups = df_all_act.columns[df_all_act.columns.duplicated()].tolist()
    if all_dups:
        logger.info(
            "activity_duplicate_columns_removed", source="all", columns=all_dups
        )
        df_all_act = df_all_act.loc[:, ~df_all_act.columns.duplicated()]

    concat = pd.concat([df_same_act, df_all_act], ignore_index=True, sort=False)
    concat = concat.drop(
        columns=list(ACTIVITY_DROP_COLS & set(concat.columns)), errors="ignore"
    )
    concat = concat.loc[:, ~concat.columns.duplicated()]
    if concat.shape[1]:
        df_activity = concat.drop_duplicates(subset="activity_chembl_id", keep="first")
    else:
        df_activity = concat
    if dict_dir is not None:
        df_activity = process_activity_table(df_activity, dict_dir, targets_path)
    logger.info(
        "activity_rows_summary",
        rows_same=len(df_same_act),
        rows_all=len(df_all_act),
        rows_concat=len(concat),
        rows_after_dedup=len(df_activity),
    )
    remaining = ACTIVITY_DROP_COLS & set(df_activity.columns)
    if remaining:
        raise ValueError(
            f"activity table still contains columns: {', '.join(sorted(remaining))}"
        )
    combined["activity"] = df_activity

    # --- pairs ------------------------------------------------------------
    df_pairs_same = normalize_pair_columns(
        add_pair_metric_columns(unify_dtypes(same["pairs_same_document"]))
    )
    df_pairs = normalize_pair_columns(
        add_pair_metric_columns(unify_dtypes(all_["pairs"]))
    )
    logger.info("pair_rows_count", table="pairs_same_document", rows=len(df_pairs_same))
    logger.info("pair_rows_count", table="pairs", rows=len(df_pairs))
    combined["pairs_same_document"] = df_pairs_same

    if "INDEPENDENT" in df_pairs.columns:
        indep_series = _safe_to_bool(df_pairs["INDEPENDENT"], "INDEPENDENT").fillna(
            False
        )
        df_pairs_independent = df_pairs[indep_series].copy()
        df_pairs_non_independent = df_pairs[~indep_series].copy()
    else:
        logger.warning("pair_independent_column_missing")
        df_pairs_independent = df_pairs.iloc[0:0].copy()
        df_pairs_non_independent = df_pairs.iloc[0:0].copy()

    combined["pairs_independent"] = df_pairs_independent
    combined["pairs_non_independent"] = df_pairs_non_independent

    return combined


def save_tables(
    tables: TableDict,
    out_dir: Path,
    cfg: Config,
    fmt: str = "csv",
) -> dict[str, Path]:
    """Persist combined tables to ``out_dir``.

    Tables with specific name suffixes are written to dedicated
    subdirectories, e.g. ``pairs_same_document`` is stored in
    ``same_document/pairs_same_document.csv``. The following suffixes are
    recognised:

    - ``*_same_document`` → ``same_document/``
    - ``*_independent`` → ``independent/``
    - ``*_non_independent`` → ``non_independent/``

    Duplicate column names are removed prior to writing to avoid ambiguous
    headers. The dropped columns are listed in a warning message to aid
    debugging.

    Parameters
    ----------
    tables:
        Mapping of entity names to :class:`pandas.DataFrame` instances.
    out_dir:
        Destination directory. Created if missing.
    cfg:
        Global configuration providing I/O settings used when writing files.
    fmt:
        Output format. Only ``"csv"`` is supported.

    Returns
    -------
    dict
        Mapping of entity names to written file paths.
    """
    if fmt != "csv":
        raise ValueError("only csv output is supported")

    paths: dict[str, Path] = {}
    for entity, df in tables.items():
        # Determine subdirectory based on table type.
        if entity.endswith("_non_independent"):
            sub_dir = out_dir / "non_independent"
        elif entity.endswith("_independent"):
            sub_dir = out_dir / "independent"
        elif entity.endswith("_same_document"):
            sub_dir = out_dir / "same_document"
        else:
            sub_dir = out_dir

        sub_dir.mkdir(parents=True, exist_ok=True)
        path = sub_dir / f"{entity}.csv"

        chembl_id_map = {"testitem": "molecule_chembl_id"}
        chembl_id_col = chembl_id_map.get(entity, f"{entity}_chembl_id")

        # Drop duplicated columns before writing to avoid ambiguous headers.
        duplicate_cols = df.columns[df.columns.duplicated()].tolist()
        if duplicate_cols:
            logger.warning(
                "duplicate_columns_removed", entity=entity, columns=duplicate_cols
            )
            df = df.loc[:, ~df.columns.duplicated()]
        # Delegate writing to the shared I/O helper to ensure consistent
        # encoding and creation of metadata sidecars.
        write_csv(df, path, cfg=cfg, col_order=[chembl_id_col])
        logger.info("file_written", rows=len(df), path=str(path))
        paths[entity] = path
    return paths


# Status processing -----------------------------------------------------------


@dataclass(frozen=True)
class StatusAPI:
    """Container for helpers related to the status configuration.

    Attributes
    ----------
    table:
        Raw dataframe sorted by ``order``.
    status_list:
        Labels ordered by ``order``.
    conditions:
        List of condition fields where ``condition_value`` is not ``"null"``.
    order_map:
        Mapping of labels to their ``order``.
    score_map:
        Mapping of labels to ``score``.

    """

    table: pd.DataFrame
    status_list: list[str]
    conditions: list[str]
    order_map: dict[str, int]
    score_map: dict[str, int]

    def pair(self, s1: str, s2: str) -> str:
        """Return the lower-ranked value between ``s1`` and ``s2``.

        Parameters
        ----------
        s1, s2:
            Values to compare.

        Returns
        -------
        str
            The label with the smaller ``order`` value.

        """
        return self.min_status([s1, s2])

    def next(self, status: str) -> str:
        """Return the label following ``status`` in ``status_list``.

        Parameters
        ----------
        status:
            Current value.

        Returns
        -------
        str
            The next label or the last element if ``status`` is unknown.

        """
        idx = (
            self.status_list.index(status)
            if status in self.status_list
            else len(self.status_list) - 1
        )
        return self.status_list[min(idx + 1, len(self.status_list) - 1)]

    def min_status(self, statuses: Iterable[str]) -> str:
        """Return the lowest ``order`` value among the provided labels.

        Parameters
        ----------
        statuses:
            Iterable of values.

        Returns
        -------
        str
            The label with the minimum ``order`` or the first element of
            ``status_list`` if none are valid.

        """
        valid = [s for s in statuses if s in self.order_map]
        if not valid:
            return self.status_list[0]
        return min(valid, key=lambda s: self.order_map[s])

    def max_status(self, statuses: Iterable[str]) -> str:
        """Return the highest ``order`` value among the provided labels.

        Parameters
        ----------
        statuses:
            Iterable of values.

        Returns
        -------
        str
            The label with the maximum ``order`` or the last element of
            ``status_list`` if none are valid.

        """
        valid = [s for s in statuses if s in self.order_map]
        if not valid:
            return self.status_list[-1]
        return max(valid, key=lambda s: self.order_map[s])

    def get_order(self, status: str) -> int:
        """Return the numeric ``order`` for the given label.

        Parameters
        ----------
        status:
            Value to resolve.

        Returns
        -------
        int
            ``order`` associated with the label; uses the last element's value
            as fallback.

        """
        return self.order_map.get(status, self.order_map[self.status_list[-1]])

    def get_score(self, status: str) -> int:
        """Return the ``score`` associated with the given label.

        Parameters
        ----------
        status:
            Value to look up.

        Returns
        -------
        int
            Score for the label; defaults to ``0`` when unknown.

        """
        return self.score_map.get(status, 0)

    def ascending(self, s1: str, s2: str) -> bool:
        """Return ``True`` if ``s1`` precedes ``s2`` by ``order`` value.

        Parameters
        ----------
        s1, s2:
            Status values to compare.

        Returns
        -------
        bool
            ``True`` when ``s1`` has a lower ``order`` than ``s2``.

        """
        return self.get_order(s1) < self.get_order(s2)

    def descending(self, s1: str, s2: str) -> bool:
        """Return ``True`` if ``s1`` follows ``s2`` by ``order`` value.

        Parameters
        ----------
        s1, s2:
            Status values to compare.

        Returns
        -------
        bool
            ``True`` when ``s1`` has a higher ``order`` than ``s2``.

        """
        return self.get_order(s1) > self.get_order(s2)

    def Next(self, status: str) -> str:
        """Compatibility alias for :meth:`next`.

        Parameters
        ----------
        status:
            Current status value.

        Returns
        -------
        str
            The next status according to :attr:`status_list`.

        """
        return self.next(status)

    def MinStatus(self, statuses: Iterable[str]) -> str:
        """Compatibility alias for :meth:`min_status`.

        Parameters
        ----------
        statuses:
            Iterable of status values.

        Returns
        -------
        str
            Result of :meth:`min_status`.

        """
        return self.min_status(statuses)

    def MaxStatus(self, statuses: Iterable[str]) -> str:
        """Compatibility alias for :meth:`max_status`.

        Parameters
        ----------
        statuses:
            Iterable of status values.

        Returns
        -------
        str
            Result of :meth:`max_status`.

        """
        return self.max_status(statuses)

    def get_min(self, statuses: Iterable[str]) -> str:
        """Return the lowest ranked status among ``statuses``.

        This method preserves backwards compatibility with code that
        previously relied on a ``get_min`` helper. When none of the provided
        values match a known status label, the first element of
        :attr:`status_list` is returned, mirroring :meth:`min_status`.

        Parameters
        ----------
        statuses:
            Iterable of status values.

        Returns
        -------
        str
            Result of :meth:`min_status`.
        """
        return self.min_status(statuses)

    def get_max(self, statuses: Iterable[str]) -> str:
        """Return the highest ranked status among ``statuses``.

        The function mirrors the behaviour of :meth:`max_status` but uses a
        name compatible with legacy code. Unknown values are ignored; if none
        of the provided statuses are recognised the last element of
        :attr:`status_list` is returned instead of raising an exception.

        Parameters
        ----------
        statuses:
            Iterable of status values.

        Returns
        -------
        str
            Result of :meth:`max_status`.
        """
        return self.max_status(statuses)

    def Ascending(self, s1: str, s2: str) -> bool:
        """Compatibility alias for :meth:`ascending`.

        Parameters
        ----------
        s1, s2:
            Status values to compare.

        Returns
        -------
        bool
            Result of :meth:`ascending`.

        """
        return self.ascending(s1, s2)

    def Descending(self, s1: str, s2: str) -> bool:
        """Compatibility alias for :meth:`descending`.

        Parameters
        ----------
        s1, s2:
            Status values to compare.

        Returns
        -------
        bool
            Result of :meth:`descending`.

        """
        return self.descending(s1, s2)


def normalize_pair_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize activity ID column names in pair tables.

    The Excel sources occasionally vary in the casing or use of underscores
    for ``activity_chembl_id1`` and ``activity_chembl_id2``. This helper
    normalises these column names so downstream processing can rely on a
    consistent schema.

    Parameters
    ----------
    df:
        DataFrame potentially containing activity ID columns with alternative
        spellings.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with recognised activity ID columns renamed to
        ``activity_id1`` and ``activity_id2``.

    """
    rename: dict[str, str] = {}
    for col in df.columns:
        key = col.lower().replace("_", "")
        if key == "activityid1":
            rename[col] = "activity_chembl_id1"
        elif key == "activityid2":
            rename[col] = "activity_chembl_id2"
    return df.rename(columns=rename)
