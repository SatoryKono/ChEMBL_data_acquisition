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


def postprocess_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and reshape merged target information.

    Parameters
    ----------
    df:
        DataFrame produced by the ``all`` pipeline in
        :mod:`get_target_data`.

    Returns
    -------
    pandas.DataFrame
        Normalised table ready for export.
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
    """

    df = pd.read_csv(input_path, sep=sep, encoding=encoding, dtype=str)
    processed = postprocess_targets(df)
    processed.to_csv(output_path, index=False, sep=sep, encoding=encoding)