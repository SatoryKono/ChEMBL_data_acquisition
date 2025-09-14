"""Schema definitions for target data."""

from __future__ import annotations

import pandera.pandas as pa

# Explicit column order for the targets table.
#
# Pandera sorts columns alphabetically when exposing ``DataFrameSchema.columns``.
# To preserve the intended ordering in data exports, the sequence is captured in
# ``TARGETS_COLUMN_ORDER`` and referenced by post-processing utilities.
TARGETS_COLUMN_ORDER: list[str] = [
    "target_chembl_id",
    "uniprotkb_Id",
    "recommended_name",
    "synonyms",
    "type",
    "uniprot_id",
    "secondary_uniprot_id",
    "gene_name",
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
]


# Definition of the schema describing the targets table used in exports.
TargetsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "target_chembl_id": pa.Column(str, required=True),
        "uniprotkb_Id": pa.Column(str, required=False),
        "recommended_name": pa.Column(str, required=False),
        "synonyms": pa.Column(str, required=False),
        "type": pa.Column(str, required=False),
        "uniprot_id": pa.Column(str, required=False),
        "secondary_uniprot_id": pa.Column(str, required=False),
        "gene_name": pa.Column(str, required=False),
        "genus": pa.Column(str, required=False),
        "superkingdom": pa.Column(str, required=False),
        "phylum": pa.Column(str, required=False),
        "taxon_id": pa.Column(int, required=False),
        "ec_number": pa.Column(str, required=False),
        "hgnc_name": pa.Column(str, required=False),
        "hgnc_id": pa.Column(str, required=False),
        "molecular_function": pa.Column(str, required=False),
        "cellular_component": pa.Column(str, required=False),
        "subcellular_location": pa.Column(str, required=False),
        "topology": pa.Column(str, required=False),
        "transmembrane": pa.Column(bool, required=False),
        "intramembrane": pa.Column(bool, required=False),
        "glycosylation": pa.Column(bool, required=False),
        "lipidation": pa.Column(bool, required=False),
        "disulfide_bond": pa.Column(bool, required=False),
        "modified_residue": pa.Column(bool, required=False),
        "phosphorylation": pa.Column(bool, required=False),
        "acetylation": pa.Column(bool, required=False),
        "ubiquitination": pa.Column(bool, required=False),
        "signal_peptide": pa.Column(bool, required=False),
        "propeptide": pa.Column(bool, required=False),
        "isoform_names": pa.Column(str, required=False),
        "isoform_ids": pa.Column(str, required=False),
        "isoform_synonyms": pa.Column(str, required=False),
        "reactions": pa.Column(str, required=False),
        "target_id": pa.Column(str, required=False),
        "IUPHAR_family_id": pa.Column(str, required=False),
        "IUPHAR_type": pa.Column(str, required=False),
        "IUPHAR_class": pa.Column(str, required=False),
        "IUPHAR_subclass": pa.Column(str, required=False),
        "IUPHAR_chain": pa.Column(str, required=False),
        "full_id_path": pa.Column(str, required=False),
        "full_name_path": pa.Column(str, required=False),
        "GuidetoPHARMACOLOGY": pa.Column(str, required=False),
    },
    ordered=True,
)

"""pa.DataFrameSchema: Validation schema for targets."""
