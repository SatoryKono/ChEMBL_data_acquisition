"""Data transformation utilities for target tables.

This module implements post-processing helpers that reshape and clean the
merged target information produced by :mod:`get_target_data`.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from library.schemas.targets import CELLULARITY_COLUMN_NAME, TARGETS_COLUMN_ORDER

from ...common.csv_utils import write_csv_deterministic
from ...common.log import logger
from ...config import Config, IoCfg
from ...postprocessing import helpers as postprocessing_helpers
from . import organism_classification

# Columns removed in the final export
REMOVE_COLUMNS: list[str] = []

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
    "protein_name_alt",
    "gtop_synonyms",
    "gene_symbol",
    "protein_synonym_list",
    "gene_symbol_list",
]

# Columns treated as text in the final table
TEXT_COLUMNS: list[str] = [
    "target_chembl_id",
    "uniprotkb_Id",
    "uniProtkbIdFallback",
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
    "protein_name_canonical",
    "protein_name_alt",
    "organism",
    "lineage_superkingdom",
    "lineage_phylum",
    "lineage_class",
    "features_topology",
    "gtop_synonyms",
    "pfam",
    "interpro",
    "xref_pdb",
    "xref_alphafold",
    "pref_name",
    CELLULARITY_COLUMN_NAME,
    "tax_id",
    "species_group_flag",
    "target_components",
    "protein_classifications",
    "cross_references",
    "gene_symbol_list",
    "protein_synonym_list",
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

TARGET_WORKBOOK_SCHEMA: dict[str, str] = {column: "Text" for column in TEXT_COLUMNS}
for column in INT_COLUMNS:
    TARGET_WORKBOOK_SCHEMA[column] = "Int64"
for column in BOOL_COLUMNS:
    TARGET_WORKBOOK_SCHEMA[column] = "Logical"

NA_MARKERS: tuple[str, ...] = ("[#N/A]",)


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


def _default_series(df: pd.DataFrame, value: str = "-") -> pd.Series:
    """Return a series of length ``len(df)`` filled with ``value``."""

    if df.empty:
        return pd.Series(dtype="string")
    return pd.Series([value] * len(df), index=df.index)


def _series_or_default(df: pd.DataFrame, column: str, default: str = "-") -> pd.Series:
    """Return ``df[column]`` with missing values replaced by ``default``."""

    if column in df.columns:
        series = df[column].copy()
        return series.fillna(default)
    return _default_series(df, default)


def _derive_genus(series: pd.Series) -> pd.Series:
    """Return the genus token extracted from an organism label series."""

    if series.empty:
        return pd.Series(dtype="string")

    as_string = series.astype("string")
    genus = as_string.str.split().str[0]
    genus = genus.where(genus.notna(), pd.NA)
    if pd.api.types.is_string_dtype(genus):
        genus = genus.where(genus.str.strip() != "", pd.NA)
    return genus.astype("string")


def _pipe_merge_columns(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    """Create a pipe-merged series from the specified ``columns``."""

    available = [col for col in columns if col in df.columns]
    if not available:
        return _default_series(df)

    merged = df.apply(
        lambda row: _pipe_merge([row.get(col) for col in available]), axis=1
    )
    merged = merged.replace("", "-")
    return merged.fillna("-")


def align_target_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Project *df* onto ``TARGETS_COLUMN_ORDER`` with derived columns."""

    aligned = pd.DataFrame(index=df.index)

    aligned["target_chembl_id"] = _series_or_default(df, "target_chembl_id")
    aligned["uniprot_id_primary"] = _series_or_default(df, "uniprotkb_Id")
    aligned["uniprot_ids_all"] = _pipe_merge_columns(
        df,
        [
            "uniprotkb_Id",
            "uniprot_id",
            "secondary_uniprot_id",
            "secondaryAccessions",
            "mapping_uniprot_id",
        ],
    )
    aligned["isoform_ids"] = _series_or_default(df, "isoform_ids")
    aligned["isoform_names"] = _series_or_default(df, "isoform_names")
    aligned["isoform_synonyms"] = _series_or_default(df, "isoform_synonyms")
    aligned["hgnc_id"] = _series_or_default(df, "hgnc_id")
    aligned["gene_symbol"] = _series_or_default(df, "gene_name")
    aligned["protein_name_canonical"] = _series_or_default(df, "recommended_name")
    aligned["protein_name_alt"] = _pipe_merge_columns(
        df, ["chembl_alternative_name", "names", "synonyms", "secondaryAccessionNames"]
    )
    aligned["organism"] = _series_or_default(df, "genus")
    aligned["taxon_id"] = _series_or_default(df, "taxon_id")
    aligned["lineage_superkingdom"] = _series_or_default(df, "superkingdom")
    aligned["lineage_phylum"] = _series_or_default(df, "phylum")
    aligned["lineage_class"] = _series_or_default(df, "lineage_class")
    aligned["sequence_length"] = _series_or_default(df, "sequence_length")
    aligned["features_signal_peptide"] = _series_or_default(df, "signal_peptide")
    aligned["features_transmembrane"] = _series_or_default(df, "transmembrane")
    aligned["features_topology"] = _series_or_default(df, "topology")
    aligned["ptm_glycosylation"] = _series_or_default(df, "glycosylation")
    aligned["ptm_lipidation"] = _series_or_default(df, "lipidation")
    aligned["ptm_disulfide_bond"] = _series_or_default(df, "disulfide_bond")
    aligned["ptm_modified_residue"] = _series_or_default(df, "modified_residue")
    aligned["xref_chembl"] = aligned["target_chembl_id"]
    aligned["xref_uniprot"] = aligned["uniprot_id_primary"]
    aligned["xref_ensembl"] = _series_or_default(df, "xref_ensembl")
    aligned["xref_iuphar"] = _series_or_default(df, "target_id")
    aligned["gtop_target_id"] = _series_or_default(df, "GuidetoPHARMACOLOGY")
    aligned["gtop_synonyms"] = _series_or_default(df, "synonyms")
    aligned["gtop_natural_ligands_n"] = _series_or_default(df, "gtop_natural_ligands_n")
    aligned["gtop_interactions_n"] = _series_or_default(df, "gtop_interactions_n")
    aligned["gtop_function_text_short"] = _series_or_default(
        df, "gtop_function_text_short"
    )
    aligned["uniprot_last_update"] = _series_or_default(df, "uniprot_last_update")
    aligned["uniprot_version"] = _series_or_default(df, "uniprot_version")
    aligned["pipeline_version"] = _series_or_default(df, "pipeline_version")
    aligned["timestamp_utc"] = _series_or_default(df, "timestamp_utc")
    aligned["pfam"] = _series_or_default(df, "Pfam")
    aligned["interpro"] = _series_or_default(df, "InterPro")
    aligned["xref_pdb"] = _series_or_default(df, "xref_pdb")
    aligned["xref_alphafold"] = _series_or_default(df, "xref_alphafold")
    aligned["hgnc_name"] = _series_or_default(df, "hgnc_name")
    aligned["uniProtkbId"] = _series_or_default(df, "uniProtkbId")
    aligned["secondaryAccessions"] = _series_or_default(df, "secondaryAccessions")
    aligned["recommendedName"] = _series_or_default(df, "recommendedName")
    aligned["geneName"] = _series_or_default(df, "geneName")
    aligned["secondaryAccessionNames"] = _series_or_default(
        df, "secondaryAccessionNames"
    )
    aligned["molecular_function"] = _series_or_default(df, "molecular_function")
    aligned["cellular_component"] = _series_or_default(df, "cellular_component")
    aligned["subcellular_location"] = _series_or_default(df, "subcellular_location")
    aligned["topology"] = _series_or_default(df, "topology")
    aligned["transmembrane"] = _series_or_default(df, "transmembrane")
    aligned["intramembrane"] = _series_or_default(df, "intramembrane")
    aligned["glycosylation"] = _series_or_default(df, "glycosylation")
    aligned["lipidation"] = _series_or_default(df, "lipidation")
    aligned["disulfide_bond"] = _series_or_default(df, "disulfide_bond")
    aligned["modified_residue"] = _series_or_default(df, "modified_residue")
    aligned["phosphorylation"] = _series_or_default(df, "phosphorylation")
    aligned["acetylation"] = _series_or_default(df, "acetylation")
    aligned["ubiquitination"] = _series_or_default(df, "ubiquitination")
    aligned["signal_peptide"] = _series_or_default(df, "signal_peptide")
    aligned["propeptide"] = _series_or_default(df, "propeptide")
    aligned["GuidetoPHARMACOLOGY"] = _series_or_default(df, "GuidetoPHARMACOLOGY")
    aligned["family"] = _series_or_default(df, "family")
    aligned["SUPFAM"] = _series_or_default(df, "SUPFAM")
    aligned["PROSITE"] = _series_or_default(df, "PROSITE")
    aligned["InterPro"] = _series_or_default(df, "InterPro")
    aligned["Pfam"] = _series_or_default(df, "Pfam")
    aligned["PRINTS"] = _series_or_default(df, "PRINTS")
    aligned["TCDB"] = _series_or_default(df, "TCDB")
    aligned["pref_name"] = _series_or_default(df, "pref_name")
    aligned[CELLULARITY_COLUMN_NAME] = _series_or_default(df, CELLULARITY_COLUMN_NAME)
    aligned["tax_id"] = _series_or_default(df, "tax_id")
    aligned["species_group_flag"] = _series_or_default(df, "species_group_flag")
    aligned["target_components"] = _series_or_default(df, "target_components")
    aligned["protein_classifications"] = _series_or_default(
        df, "protein_classifications"
    )
    aligned["cross_references"] = _series_or_default(df, "cross_references")
    aligned["gene_symbol_list"] = _series_or_default(df, "gene")
    aligned["protein_synonym_list"] = _series_or_default(df, "synonyms")
    aligned["reactions"] = _series_or_default(df, "reactions")
    aligned["reaction_ec_numbers"] = _series_or_default(df, "reaction_ec_numbers")
    aligned["protein_class_pred_L1"] = _series_or_default(df, "protein_class_pred_L1")
    aligned["protein_class_pred_L2"] = _series_or_default(df, "protein_class_pred_L2")
    aligned["protein_class_pred_L3"] = _series_or_default(df, "protein_class_pred_L3")
    aligned["protein_class_pred_rule_id"] = _series_or_default(
        df, "protein_class_pred_rule_id"
    )
    aligned["protein_class_pred_evidence"] = _series_or_default(
        df, "protein_class_pred_evidence"
    )
    aligned["protein_class_pred_confidence"] = _series_or_default(
        df, "protein_class_pred_confidence"
    )
    aligned["iuphar_target_id"] = _series_or_default(df, "target_id")
    aligned["iuphar_family_id"] = _series_or_default(df, "IUPHAR_family_id")
    aligned["iuphar_type"] = _series_or_default(df, "IUPHAR_type")
    aligned["iuphar_class"] = _series_or_default(df, "IUPHAR_class")
    aligned["iuphar_subclass"] = _series_or_default(df, "IUPHAR_subclass")
    aligned["iuphar_chain"] = _series_or_default(df, "IUPHAR_chain")
    aligned["iuphar_name"] = _series_or_default(df, "iuphar_name")
    aligned["iuphar_full_id_path"] = _series_or_default(df, "full_id_path")
    aligned["iuphar_full_name_path"] = _series_or_default(df, "full_name_path")

    aligned = aligned.reindex(columns=TARGETS_COLUMN_ORDER, fill_value="-")
    aligned = aligned.fillna("-").astype("string")
    lowercase_after_align = [
        "features_signal_peptide",
        "features_transmembrane",
        "ptm_glycosylation",
        "ptm_lipidation",
        "ptm_disulfide_bond",
        "ptm_modified_residue",
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
    for col in lowercase_after_align:
        if col in aligned.columns:
            aligned[col] = aligned[col].str.lower()
    aligned = aligned.mask(aligned == "", "-")
    return aligned


def postprocess_targets(
    df: pd.DataFrame, *, chembl_col: str = "target_chembl_id"
) -> pd.DataFrame:
    """Clean and reshape merged target information.

    The function normalises UniProt identifiers, resolves gene names and
    synonyms, fills optional columns and standardises the output ordering so
    that the resulting table is ready for downstream export.

    Parameters
    ----------
    df:
        DataFrame produced by the ``all`` pipeline in
        :mod:`get_target_data`.
    chembl_col:
        Name of the column containing ChEMBL target identifiers. The input
        column will be preserved in the returned DataFrame. Defaults to
        ``"target_chembl_id"``.

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
    internal_id = "target_chembl_id"
    if chembl_col in df.columns and chembl_col != internal_id:
        df = df.rename(columns={chembl_col: internal_id})

    # --- normalise identifiers -------------------------------------------------
    raw_uniprot_ids = df.get("uniProtkbId")
    if raw_uniprot_ids is not None:
        raw_uniprot_ids = raw_uniprot_ids.astype("string")
    else:
        raw_uniprot_ids = pd.Series(pd.NA, index=df.index, dtype="string")

    fallback_uniprot_ids = df.get("uniProtkbIdFallback")
    if fallback_uniprot_ids is not None:
        fallback_uniprot_ids = fallback_uniprot_ids.astype("string")
        resolved_uniprot_ids = raw_uniprot_ids.fillna(fallback_uniprot_ids)
    else:
        resolved_uniprot_ids = raw_uniprot_ids.copy()

    uniprot_id_values = df.get("uniprot_id")
    if uniprot_id_values is not None:
        resolved_uniprot_ids = resolved_uniprot_ids.fillna(
            uniprot_id_values.astype("string")
        )

    df["uniprotkb_Id"] = (
        resolved_uniprot_ids.fillna("")
        .astype(str)
        .str.split("_")
        .str[0]
        .str.split("-")
        .str[0]
    )
    secondary_accessions = df.get("secondaryAccessions")
    if secondary_accessions is not None:
        secondary_accessions = secondary_accessions.astype("string")
    else:
        secondary_accessions = pd.Series(pd.NA, index=df.index, dtype="string")

    if "uniprot_id" in df.columns:
        secondary_fallback = df["uniprot_id"].astype("string")
    else:
        secondary_fallback = pd.Series(pd.NA, index=df.index, dtype="string")

    df["secondary_uniprot_id"] = secondary_accessions.fillna(secondary_fallback)

    # Copy "pref_name" into the normalised "recommended_name" column but keep the
    # original field for downstream exports
    if "pref_name" in df.columns and "recommended_name" not in df.columns:
        df["recommended_name"] = df["pref_name"]

    # --- gene name handling -----------------------------------------------------
    df["gene_name_x"] = df.get("gene_name_x", pd.Series(dtype=str)).replace(
        {"51.1rMVA_034": "N1L"}
    )
    df["gene_name"] = df.get("geneName", pd.Series(dtype=str))

    if "gene" in df.columns:
        df["gene"] = df["gene"].astype("string").fillna("")
    else:
        df["gene"] = pd.Series("", index=df.index, dtype="string")
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
    optional_columns = [
        "isoform_names",
        "isoform_ids",
        "isoform_synonyms",
        "reactions",
        "full_id_path",
        "full_name_path",
    ]
    existing_optional = [col for col in optional_columns if col in df.columns]
    if existing_optional:
        df[existing_optional] = df[existing_optional].fillna("-")

    missing_optional = [col for col in optional_columns if col not in df.columns]
    if missing_optional:
        df[missing_optional] = "-"

    # --- deduplicate gene field -------------------------------------------------
    df["gene"] = df["gene"].apply(lambda v: _pipe_merge([v]))

    # --- final column ordering --------------------------------------------------
    schema_cols = [
        internal_id if col == "target_chembl_id" else col
        for col in TARGETS_COLUMN_ORDER
    ]
    extra_cols = sorted(c for c in df.columns if c not in schema_cols)
    ordered_cols = schema_cols + extra_cols
    for col in ordered_cols:
        if col not in df.columns:
            df[col] = "-"
    return df[ordered_cols].rename(columns={internal_id: chembl_col})


def postprocess_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    cfg: IoCfg,
    chembl_col: str = "target_chembl_id",
    sep: str | None = None,
    encoding: str | None = None,
) -> None:
    """Read a CSV, post-process and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by ``scripts/get_target_data.py all``.
    output_path:
        Destination path for the cleaned CSV file.
    cfg:
        I/O configuration providing default CSV parameters.
    chembl_col:
        Name of the column containing ChEMBL target identifiers. Passed to
        :func:`postprocess_targets`.
    sep:
        Field delimiter of the CSV files. Defaults to ``cfg.csv_sep``.
    encoding:
        Text encoding of the CSV files. Defaults to ``cfg.csv_encoding``.

    Tests
    -----
    Verified in
    :func:`tests.test_target_postprocessing.test_postprocess_file_roundtrip`.

    """
    sep = sep or cfg.csv_sep
    encoding = encoding or cfg.csv_encoding
    schema = dict(TARGET_WORKBOOK_SCHEMA)
    schema.setdefault(chembl_col, "Text")

    encodings: list[str] = []
    if encoding:
        encodings.append(encoding)
    encodings.extend(postprocessing_helpers.ENCODING_FALLBACKS)

    df = postprocessing_helpers.read_csv_strict(
        input_path,
        encoding=encodings,
        dtype_map=schema,
        na_values=NA_MARKERS,
        separators=(sep,),
    )
    processed = postprocess_targets(df, chembl_col=chembl_col)
    write_csv_deterministic(
        processed,
        output_path,
        col_order=list(processed.columns),
        key_cols=[chembl_col],
        sep=sep,
        encoding=encoding,
        cfg=None,
    )


def finalise_targets(
    df: pd.DataFrame,
    *,
    chembl_col: str = "target_chembl_id",
    uniprot_col: str = "uniprotkb_Id",
    genus_col: str = "genus",
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str = "lineage_class",
) -> pd.DataFrame:
    """Apply final cleaning steps and derive organism cellularity.

    The function mirrors the original Power Query pipeline but now infers the
    organism ``type`` directly from UniProt taxonomy columns using
    :func:`organism_classification.add_cellularity_smart`. It filters rows with
    missing identifiers, enforces data types and normalises selected text
    fields before aligning the output schema.

    Parameters
    ----------
    df:
        DataFrame produced by :func:`postprocess_targets`.
    chembl_col:
        Column containing ChEMBL target identifiers.
    uniprot_col:
        Column containing UniProt identifiers.
    genus_col:
        Column containing the organism genus used for cellularity inference.
    superkingdom_col, phylum_col, class_col:
        Taxonomy lineage columns required by the cellularity classifier.

    Returns
    -------
    pandas.DataFrame
        Cleaned table ready for export.

    Notes
    -----
    If ``df`` already contains a ``type`` column, it will be renamed to
    ``AddCellularitySmart `` (the original Power Query export name) before the
    new classification is appended to avoid conflicting suffixes.
    """

    df = df.copy()
    internal_mapping = {
        chembl_col: "target_chembl_id",
        uniprot_col: "uniprotkb_Id",
        genus_col: "genus",
        superkingdom_col: "lineage_superkingdom",
        phylum_col: "lineage_phylum",
        class_col: "lineage_class",
    }
    df = df.rename(columns={k: v for k, v in internal_mapping.items() if k != v})

    genus_missing_mask: pd.Series
    if "genus" in df.columns:
        genus_series = df["genus"].astype("string")
        genus_missing_mask = genus_series.isna() | genus_series.str.strip().isin(
            ["", "nan"]
        )
    else:
        genus_missing_mask = pd.Series(True, index=df.index)

    if genus_missing_mask.any():
        fallback_candidates = [
            col
            for col in (
                "organism",
                "organism_name",
                "organism_x",
                "organism_y",
            )
            if col in df.columns
        ]
        for candidate in fallback_candidates:
            derived = _derive_genus(df[candidate])
            if "genus" in df.columns:
                to_fill = genus_missing_mask & derived.notna()
                if to_fill.any():
                    df.loc[to_fill, "genus"] = derived[to_fill]
            else:
                df["genus"] = derived

            genus_series = df["genus"].astype("string")
            genus_missing_mask = genus_series.isna() | genus_series.str.strip().isin(
                ["", "nan"]
            )
            if not genus_missing_mask.any():
                logger.debug(
                    "Filled missing 'genus' values from '%s' column", candidate
                )
                break

    if "lineage_superkingdom" not in df.columns and "superkingdom" in df.columns:
        df["lineage_superkingdom"] = df["superkingdom"]
    if "lineage_phylum" not in df.columns and "phylum" in df.columns:
        df["lineage_phylum"] = df["phylum"]
    if "lineage_class" not in df.columns and "class" in df.columns:
        df["lineage_class"] = df["class"]

    required_columns = [
        "target_chembl_id",
        "uniprotkb_Id",
        "genus",
        "lineage_superkingdom",
        "lineage_phylum",
        "lineage_class",
    ]
    _validate_columns(df, required_columns)

    uniprot_series = df["uniprotkb_Id"].astype("string")
    fallback_series = df.get("uniProtkbIdFallback")
    if fallback_series is not None:
        fallback_series = fallback_series.astype("string")
    else:
        fallback_series = pd.Series(pd.NA, index=df.index, dtype="string")

    missing_uniprot = (
        uniprot_series.isna()
        | (uniprot_series.str.strip() == "")
        | (uniprot_series.str.lower() == "nan")
    )
    fallback_available = ~(
        fallback_series.isna()
        | (fallback_series.str.strip() == "")
        | (fallback_series.str.lower() == "nan")
    )

    mask_nan = missing_uniprot & ~fallback_available
    if mask_nan.any():
        logger.debug("Dropping %d rows with missing UniProt IDs", mask_nan.sum())
    df = df[~mask_nan]

    before = len(df)
    df = df.drop_duplicates(subset="target_chembl_id", keep="first")
    logger.debug("Removed %d duplicate %s rows", before - len(df), chembl_col)

    for col in TEXT_COLUMNS:
        if col in df.columns:
            df[col] = df[col].astype("string")
    # for col in INT_COLUMNS:
    #     if col in df.columns:
    #         df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    # for col in BOOL_COLUMNS:
    #     if col in df.columns:
    #         df[col] = (
    #             df[col]
    #             .astype("string")
    #             .str.lower()
    #             .map({"true": True, "false": False})
    #             .astype("boolean")
    #         )

    if "type" in df.columns:
        logger.debug("Renaming existing 'type' column to '%s'", CELLULARITY_COLUMN_NAME)
        df = df.rename(columns={"type": CELLULARITY_COLUMN_NAME})

    df = organism_classification.add_cellularity_smart(
        df,
        genus_col="genus",
        superkingdom_col="lineage_superkingdom",
        phylum_col="lineage_phylum",
        class_col="lineage_class",
        output_col=CELLULARITY_COLUMN_NAME,
    )

    df = df.drop(columns=[c for c in REMOVE_COLUMNS if c in df.columns])

    for col in LOWERCASE_COLUMNS:
        if col in df.columns:
            df[col] = df[col].astype("string").str.lower()

    df = align_target_columns(df)
    return df


def finalise_file(
    input_path: Path | str,
    output_path: Path | str,
    *,
    cfg: Config,
    chembl_col: str = "target_chembl_id",
    uniprot_col: str = "uniprotkb_Id",
    genus_col: str = "genus",
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str = "lineage_class",
    sep: str | None = None,
    encoding: str | None = None,
) -> None:
    """Read CSV files, finalise the target table and write the result.

    Parameters
    ----------
    input_path:
        Path to the CSV file produced by :func:`postprocess_file`.
    output_path:
        Destination path for the cleaned CSV file.
    cfg:
        Application configuration providing CSV defaults and resource paths.
    chembl_col:
        Column containing ChEMBL target identifiers.
    uniprot_col:
        Column containing UniProt identifiers.
    genus_col:
        Column containing the organism genus used for merging.
    superkingdom_col, phylum_col, class_col:
        Taxonomy lineage columns required by the cellularity classifier.
    sep:
        Field delimiter of the CSV files. Defaults to ``cfg.io.csv_sep``.
    encoding:
        Text encoding of the CSV files. Defaults to ``cfg.io.csv_encoding``.

    """
    sep = sep or cfg.io.csv_sep
    encoding = encoding or cfg.io.csv_encoding
    schema = dict(TARGET_WORKBOOK_SCHEMA)
    schema.setdefault(chembl_col, "Text")

    encodings: list[str] = []
    if encoding:
        encodings.append(encoding)
    encodings.extend(postprocessing_helpers.ENCODING_FALLBACKS)

    df = postprocessing_helpers.read_csv_strict(
        input_path,
        encoding=encodings,
        dtype_map=schema,
        na_values=NA_MARKERS,
        separators=(sep,),
    )
    processed = finalise_targets(
        df,
        chembl_col=chembl_col,
        uniprot_col=uniprot_col,
        genus_col=genus_col,
        superkingdom_col=superkingdom_col,
        phylum_col=phylum_col,
        class_col=class_col,
    )
    write_csv_deterministic(
        processed,
        output_path,
        col_order=list(processed.columns),
        key_cols=[chembl_col],
        sep=sep,
        encoding=encoding,
        cfg=cfg,
    )
