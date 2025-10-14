"""Helpers for deriving the multifunctional enzyme flag."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import pandas as pd

_COLUMNS_TO_REMOVE: tuple[str, ...] = (
    "isoform_ids",
    "isoform_names",
    "isoform_synonyms",
    "organism",
    "taxon_id",
    "lineage_superkingdom",
    "lineage_phylum",
    "lineage_class",
    "xref_chembl",
    "gtop_synonyms",
    "gtop_natural_ligands_n",
    "gtop_interactions_n",
    "gtop_function_text_short",
    "uniProtkbId",
    "hgnc_name",
    "secondaryAccessionNames",
    "transmembrane",
    "intramembrane",
    "hgnc_id",
    "features_signal_peptide",
    "tax_id",
    "species_group_flag",
    "family",
    "xref_uniprot",
    "xref_ensembl",
    "pfam",
    "interpro",
    "xref_pdb",
    "xref_alphafold",
    "SUPFAM",
    "PROSITE",
    "InterPro",
    "Pfam",
    "PRINTS",
    "TCDB",
    "target_type",
    "protein_classifications",
    "protein_name_alt",
    "recommendedName",
    "pref_name",
    "target_components",
    "cross_references",
    "gene_symbol_list",
    "protein_synonym_list",
    "ptm_glycosylation",
    "ptm_lipidation",
    "ptm_disulfide_bond",
    "ptm_modified_residue",
    "glycosylation",
    "lipidation",
    "disulfide_bond",
    "modified_residue",
    "phosphorylation",
    "acetylation",
    "ubiquitination",
    "signal_peptide",
    "propeptide",
    "gene_symbol",
    "protein_name_canonical",
    "sequence_length",
    "features_transmembrane",
    "features_topology",
    "xref_iuphar",
    "gtop_target_id",
    "uniprot_last_update",
    "uniprot_version",
    "pipeline_version",
    "timestamp_utc",
    "secondaryAccessions",
    "geneName",
    "molecular_function",
    "cellular_component",
    "subcellular_location",
    "topology",
    "GuidetoPHARMACOLOGY",
    "protein_class_pred_L1",
    "protein_class_pred_L2",
    "protein_class_pred_L3",
    "protein_class_pred_rule_id",
    "protein_class_pred_evidence",
    "protein_class_pred_confidence",
    "iuphar_target_id",
    "iuphar_family_id",
    "iuphar_type",
    "iuphar_class",
    "iuphar_subclass",
    "iuphar_chain",
    "iuphar_name",
    "iuphar_full_id_path",
    "iuphar_full_name_path",
)


def _distinct(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def _transform_reaction_ec_numbers(value: Any) -> list[str]:
    if value is None or value is pd.NA:
        return []
    try:
        if isinstance(value, float) and pd.isna(value):
            return []
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass
    text = str(value)
    parts = text.split("|")
    parts = _distinct(parts)
    prefixes: list[str] = []
    for part in parts:
        prefix = part.split(".")[0] if part else ""
        prefixes.append(prefix)
    prefixes = _distinct(prefixes)
    return prefixes


def compute_multifunctional(source: pd.DataFrame) -> pd.DataFrame:
    """Replicate the ``multifunctional`` helper from the M script."""

    # Some datasets omit metadata fields such as ``target_type``.  We therefore
    # trim only the columns that actually exist instead of relying on
    # ``errors="ignore"`` which has proved brittle across pandas versions.
    removable_columns = [
        column for column in _COLUMNS_TO_REMOVE if column in source.columns
    ]
    if removable_columns:
        trimmed = source.drop(columns=removable_columns)
    else:
        trimmed = source.copy()
    result = trimmed.copy()

    if "target_chembl_id" not in result.columns:
        fill_values = [""] * len(result.index)
        result.insert(
            0,
            "target_chembl_id",
            pd.Series(fill_values, index=result.index, dtype="string"),
        )
    else:
        result["target_chembl_id"] = result["target_chembl_id"].astype("string")

    default_reaction_numbers = [pd.NA] * len(result.index)
    reaction_numbers = result.get(
        "reaction_ec_numbers",
        pd.Series(default_reaction_numbers, index=result.index, dtype="string"),
    )
    result["reaction_ec_numbers"] = reaction_numbers
    result["reaction_ec_numbers"] = result["reaction_ec_numbers"].map(
        _transform_reaction_ec_numbers
    )
    result["multifunctional_enzyme"] = result["reaction_ec_numbers"].map(
        lambda values: len(values) > 1
    )
    return result


__all__ = ["compute_multifunctional"]
