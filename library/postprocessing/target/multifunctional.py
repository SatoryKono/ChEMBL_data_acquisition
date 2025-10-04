"""Multifunctional enzyme detection logic translated from Power Query."""

from __future__ import annotations

from typing import Any

import pandas as pd

from .helpers import distinct_preserve_order

__all__ = ["compute_multifunctional_table"]

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


def _transform_ec_numbers(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value)
    if text == "":
        return []
    segments = text.split("|")
    distinct_segments = distinct_preserve_order(segments)
    prefixes = [segment.split(".")[0] if segment else "" for segment in distinct_segments]
    distinct_prefixes = distinct_preserve_order(prefixes)
    return distinct_prefixes


def compute_multifunctional_table(source: pd.DataFrame) -> pd.DataFrame:
    """Replicate the Power Query ``multifunctional`` helper."""

    reduced = source.drop(columns=_COLUMNS_TO_REMOVE, errors="ignore").copy()
    if "reaction_ec_numbers" not in reduced.columns:
        reduced["reaction_ec_numbers"] = [
            [] for _ in range(len(reduced))
        ]
    reduced["reaction_ec_numbers"] = reduced["reaction_ec_numbers"].apply(_transform_ec_numbers)
    reduced["multifunctional_enzyme"] = reduced["reaction_ec_numbers"].apply(
        lambda values: len(values) > 1
    )
    return reduced
