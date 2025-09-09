"""Data transformation utilities for target tables.

This module implements post-processing helpers that reshape and clean the
merged target information produced by :mod:`get_target_data`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import logging

import pandas as pd

logger = logging.getLogger(__name__)


# Columns removed in the final export
REMOVE_COLUMNS: list[str] = [
    "SUPFAM",
    "PROSITE",
    "InterPro",
    "Pfam",
    "PRINTS",
    "TCDB",
]

# Columns that should be transformed to lowercase
LOWERCASE_COLUMNS: list[str] = [
    "isoform_synonyms",
    "isoform_names",
    "topology",
    "cellular_component",
    "subcellular_location",
    "molecular_function",
    "recommended_name",
    "synonyms",
]

# Columns treated as text in the final table
TEXT_COLUMNS: list[str] = [
    "chembl_id",
    "uniprotkb_Id",
    "uniprot_id",
    "secondary_uniprot_id",
    "gene_name",
    "recommended_name",
    "synonyms",
    "genus",
    "superkingdom",
    "phylum",
    "ec_number",
    "hgnc_name",
    "hgnc_id",
    "molecular_function",
    "cellular_component",
    "subcellular_location",
    "topology",
    "isoform_names",
    "isoform_ids",
    "isoform_synonyms",
    "reactions",
    "target_id",
    "IUPHAR_family_id",
    "IUPHAR_type",
    "IUPHAR_class",
    "IUPHAR_subclass",
    "IUPHAR_chain",
    "full_id_path",
    "full_name_path",
    "GuidetoPHARMACOLOGY",
]

# Integer and boolean columns
INT_COLUMNS: list[str] = ["taxon_id"]

BOOL_COLUMNS: list[str] = [
    "transmembrane",
    "intramembrane",
    "glycosylation",
    "lipidation",
    "disulfide_bond",
    "modified_residue",
    "phosphorylation",
    "acetylation",
    "ubiquitination",
    "signal_peptide",
    "propeptide",
]


def _pipe_merge(values: Iterable[str | float | None]) -> str:
    """Return a ``"|"``-separated string of unique tokens.

    Parameters
    ----------
    values:
        Iterable of pipe-delimited strings. ``None`` and empty strings are
        ignored.

    Returns
    -------
    str
        Unique tokens joined by ``"|"`` in their first appearance order.

    """
    tokens: list[str] = []
    seen: set[str] = set()
    for value in values:
        if isinstance(value, str) and value:
            for part in (p.strip() for p in value.split("|") if p.strip()):
                if part not in seen:
                    seen.add(part)
                    tokens.append(part)
    return "|".join(tokens)


def _first_token(value: str | float | None) -> str:
    """Return the first token from a pipe-delimited string."""
    if isinstance(value, str) and value:
        return value.split("|")[0]
    return ""


def _validate_columns(df: pd.DataFrame, required: Iterable[str]) -> None:
    """Ensure that *df* contains all *required* columns."""
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"missing required columns: {', '.join(sorted(missing))}")


def postprocess_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and reshape merged target information.

    The function normalises UniProt identifiers, resolves gene names and
    synonyms, fills optional columns and standardises the output ordering so
    that the resulting table is ready for downstream export.

    Parameters
    ----------
    df:
        DataFrame produced by the ``all`` pipeline in
        :mod:`get_target_data`.

    Returns
    -------
    pandas.DataFrame
        Normalised table ready for export.

    Tests
    -----
    Behaviour is exercised in
    :func:`tests.test_target_postprocessing.test_postprocess_targets_merges_and_normalises`.

    """
    df = df.copy()

    # --- normalise identifiers -------------------------------------------------
    df["uniprotkb_Id"] = (
        df.get("uniProtkbId", pd.Series(dtype=str))
        .astype(str)
        .str.split("_")
        .str[0]
        .str.split("-")
        .str[0]
    )
    df["secondary_uniprot_id"] = df.get(
        "secondaryAccessions", pd.Series(dtype=str)
    ).fillna(df.get("uniprot_id"))

    # Rename "pref_name" to the exported "recommended_name" column
    if "pref_name" in df.columns:
        df = df.rename(columns={"pref_name": "recommended_name"})

    # --- gene name handling -----------------------------------------------------
    df["gene_name_x"] = df.get("gene_name_x", pd.Series(dtype=str)).replace(
        {"51.1rMVA_034": "N1L"}
    )
    df["gene_name"] = df.get("geneName", pd.Series(dtype=str))
    mask = df["gene_name"].isna() | (df["gene_name"] == "")
    df.loc[mask, "gene_name"] = df.loc[mask, "gene_name_x"].fillna("")
    mask = df["gene_name"] == ""
    df.loc[mask, "gene_name"] = df.loc[mask, "gene"].apply(_first_token)
    df["gene_name"] = df["gene_name"].replace("", "-").str.upper()

    df["gene"] = df.apply(
        lambda r: _pipe_merge([r.get("gene"), r.get("gene_name")]), axis=1
    )
    df["gene"] = df["gene"].replace("", "-").str.upper()

    # --- synonyms --------------------------------------------------------------
    synonym_fields = [
        "gene",
        "secondaryAccessionNames",
        "component_description",
        "chembl_alternative_name",
        "recommendedName",
        "names",
        "gene_name_x",
        "synonyms_x",
        "synonyms",
    ]
    df["synonyms"] = df.apply(
        lambda r: _pipe_merge([r.get(f) for f in synonym_fields]), axis=1
    )
    df["synonyms"] = (
        df["synonyms"]
        .str.replace("||", "|", regex=False)
        .str.replace("| ", "|", regex=False)
        .str.replace(" |", "|", regex=False)
        .str.lower()
    )

    # --- EC numbers ------------------------------------------------------------
    df["ec_number"] = df.apply(
        lambda r: _pipe_merge([r.get("ec_number"), r.get("ec_code")]), axis=1
    )
    df["ec_number"] = df["ec_number"].replace("", "-")

    # --- fill optional columns --------------------------------------------------
    for col in [
        "isoform_names",
        "isoform_ids",
        "isoform_synonyms",
        "reactions",
        "full_id_path",
        "full_name_path",
    ]:
        if col in df.columns:
            df[col] = df[col].fillna("-")
        else:
            df[col] = "-"

    # --- deduplicate gene field -------------------------------------------------
    df["gene"] = df["gene"].apply(lambda v: _pipe_merge([v]))

    # --- final column ordering --------------------------------------------------
    columns = [
        "chembl_id",
        "uniprotkb_Id",
        "uniprot_id",
        "secondary_uniprot_id",
        "gene_name",
        "recommended_name",
        "synonyms",
        "genus",
        "superkingdom",
        "phylum",
        "taxon_id",
        "ec_number",
        "hgnc_name",
        "hgnc_id",
        "molecular_function",
        "cellular_component",
        "subcellular_location",
        "topology",
        "transmembrane",
        "intramembrane",
        "glycosylation",
        "lipidation",
        "disulfide_bond",
        "modified_residue",
        "phosphorylation",
        "acetylation",
        "ubiquitination",
        "signal_peptide",
        "propeptide",
        "isoform_names",
        "isoform_ids",
        "isoform_synonyms",
        "reactions",
        "target_id",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "full_id_path",
        "full_name_path",
        "GuidetoPHARMACOLOGY",
        "SUPFAM",
        "PROSITE",
        "InterPro",
        "Pfam",
        "PRINTS",
        "TCDB",
    ]
    for col in columns:
        if col not in df.columns:
            df[col] = "-"
    return df[columns]


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> None:
    """Read a CSV, post-process and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by ``get_target_data.py all``.
    output_path:
        Destination path for the cleaned CSV file.
    sep:
        Field delimiter of the CSV files.
    encoding:
        Text encoding of the CSV files.

    Tests
    -----
    Verified in
    :func:`tests.test_target_postprocessing.test_postprocess_file_roundtrip`.

    """
    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_targets(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)


def finalise_targets(df: pd.DataFrame, organism: pd.DataFrame) -> pd.DataFrame:
    """Apply final cleaning steps and organism merge.

    This function mirrors a Power Query script that was previously used to
    prepare the exportable target table. It filters rows with missing
    identifiers, enforces data types, joins the organism classification and
    normalises selected text fields.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`postprocess_targets`.
    organism:
        Lookup table with at least ``genus`` and ``type`` columns.

    Returns
    -------
    pandas.DataFrame
        Cleaned table ready for export.

    """
    _validate_columns(df, ["chembl_id", "uniprotkb_Id", "genus"])
    _validate_columns(organism, ["genus", "type"])

    df = df.copy()

    # Drop rows where uniprotkb_Id is the string "nan"
    mask_nan = df["uniprotkb_Id"].astype(str) == "nan"
    if mask_nan.any():
        logger.debug("Dropping %d rows with missing UniProt IDs", mask_nan.sum())
    df = df[~mask_nan]

    # Remove duplicate chembl_id entries
    before = len(df)
    df = df.drop_duplicates(subset="chembl_id", keep="first")
    logger.debug("Removed %d duplicate chembl_id rows", before - len(df))

    # Enforce column types
    for col in TEXT_COLUMNS:
        if col in df.columns:
            df[col] = df[col].astype("string")
    for col in INT_COLUMNS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    for col in BOOL_COLUMNS:
        if col in df.columns:
            df[col] = (
                df[col]
                .astype("string")  # normalise mixed inputs
                .str.lower()
                .map({"true": True, "false": False})
                .astype("boolean")
            )

    # Merge organism classification and add type column
    df = df.merge(organism[["genus", "type"]], on="genus", how="left")
    df["type"] = df["type"].astype("string")

    # Remove unwanted columns
    df = df.drop(columns=[c for c in REMOVE_COLUMNS if c in df.columns])

    # Lowercase selected text columns
    for col in LOWERCASE_COLUMNS:
        if col in df.columns:
            df[col] = df[col].astype("string").str.lower()

    return df


def finalise_file(
    input_path: Path | str,
    organism_path: Path | str,
    output_path: Path | str,
    *,
    sep: str = ",",
    encoding: str = "utf8",
) -> None:
    """Read CSV files, finalise the target table and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by :func:`postprocess_file`.
    organism_path:
        Path to a CSV containing organism ``genus`` and ``type`` columns.
    output_path:
        Destination path for the cleaned CSV file.
    sep:
        Field delimiter of the CSV files.
    encoding:
        Text encoding of the CSV files.

    """
    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    organism = pd.read_csv(organism_path, sep=sep, encoding=encoding, dtype=str)
    processed = finalise_targets(df, organism)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)
