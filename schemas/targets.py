"""Schema definitions for target data."""

from __future__ import annotations

from typing import cast

import pandera.pandas as pa
from pandera.dtypes import DataType

PA_ANY = cast(DataType, None)

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
        "target_chembl_id": pa.Column(str, required=True, nullable=True),
        "uniprotkb_Id": pa.Column(str, required=False, nullable=True),
        "recommended_name": pa.Column(str, required=False, nullable=True),
        "synonyms": pa.Column(str, required=False, nullable=True),
        "type": pa.Column(str, required=False, nullable=True),
        "uniprot_id": pa.Column(str, required=False, nullable=True),
        "secondary_uniprot_id": pa.Column(str, required=False, nullable=True),
        "gene_name": pa.Column(str, required=False, nullable=True),
        "genus": pa.Column(str, required=False, nullable=True),
        "superkingdom": pa.Column(str, required=False, nullable=True),
        "phylum": pa.Column(str, required=False, nullable=True),
        "taxon_id": pa.Column(PA_ANY, required=False, nullable=True),
        "ec_number": pa.Column(str, required=False, nullable=True),
        "hgnc_name": pa.Column(str, required=False, nullable=True),
        "hgnc_id": pa.Column(str, required=False, nullable=True),
        "molecular_function": pa.Column(str, required=False, nullable=True),
        "cellular_component": pa.Column(str, required=False, nullable=True),
        "subcellular_location": pa.Column(str, required=False, nullable=True),
        "topology": pa.Column(str, required=False, nullable=True),
        "transmembrane": pa.Column(PA_ANY, required=False, nullable=True),
        "intramembrane": pa.Column(PA_ANY, required=False, nullable=True),
        "glycosylation": pa.Column(PA_ANY, required=False, nullable=True),
        "lipidation": pa.Column(PA_ANY, required=False, nullable=True),
        "disulfide_bond": pa.Column(PA_ANY, required=False, nullable=True),
        "modified_residue": pa.Column(PA_ANY, required=False, nullable=True),
        "phosphorylation": pa.Column(PA_ANY, required=False, nullable=True),
        "acetylation": pa.Column(PA_ANY, required=False, nullable=True),
        "ubiquitination": pa.Column(PA_ANY, required=False, nullable=True),
        "signal_peptide": pa.Column(PA_ANY, required=False, nullable=True),
        "propeptide": pa.Column(PA_ANY, required=False, nullable=True),
        "isoform_names": pa.Column(str, required=False, nullable=True),
        "isoform_ids": pa.Column(str, required=False, nullable=True),
        "isoform_synonyms": pa.Column(str, required=False, nullable=True),
        "reactions": pa.Column(str, required=False, nullable=True),
        "target_id": pa.Column(str, required=False, nullable=True),
        "IUPHAR_family_id": pa.Column(str, required=False, nullable=True),
        "IUPHAR_type": pa.Column(str, required=False, nullable=True),
        "IUPHAR_class": pa.Column(str, required=False, nullable=True),
        "IUPHAR_subclass": pa.Column(str, required=False, nullable=True),
        "IUPHAR_chain": pa.Column(str, required=False, nullable=True),
        "full_id_path": pa.Column(str, required=False, nullable=True),
        "full_name_path": pa.Column(str, required=False, nullable=True),
        "GuidetoPHARMACOLOGY": pa.Column(str, required=False, nullable=True),
    },
    ordered=True,
)

"""pa.DataFrameSchema: Validation schema for targets."""
