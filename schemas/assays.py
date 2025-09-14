"""Schema definitions for assay data."""

from __future__ import annotations

import pandera.pandas as pa

AssaysSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "assay_chembl_id": pa.Column(str, required=True),
        "document_chembl_id": pa.Column(str, required=True),
        "target_chembl_id": pa.Column(str, required=True),
        "aidx": pa.Column(str, required=False),
        "assay_classifications": pa.Column(str, required=False),
        "assay_group": pa.Column(str, required=False),
        "assay_category": pa.Column(str, required=False),
        "assay_organism": pa.Column(str, required=False),
        "assay_strain": pa.Column(str, required=False),
        "assay_subcellular_fraction": pa.Column(str, required=False),
        "assay_tax_id": pa.Column(str, required=False),
        "assay_tissue": pa.Column(str, required=False),
        "assay_type": pa.Column(str, required=False),
        "assay_type_description": pa.Column(str, required=False),
        "assay_test_type": pa.Column(str, required=False),
        "bao_format": pa.Column(str, required=False),
        "bao_label": pa.Column(str, required=False),
        "cell_chembl_id": pa.Column(str, required=False),
        "confidence_score": pa.Column(str, required=False),
        "description": pa.Column(str, required=False),
        "relationship_type": pa.Column(str, required=False),
        "src_assay_id": pa.Column(str, required=False),
        "src_id": pa.Column(str, required=False),
        "assay_parameters": pa.Column(str, required=False),
        "assay_cell_type": pa.Column(str, required=False),
        "tissue_chembl_id": pa.Column(str, required=False),
        "variant_sequence.isoform": pa.Column(str, required=False),
        "variant_sequence.mutation": pa.Column(str, required=False),
        "variant_sequence.sequence": pa.Column(str, required=False),
        "assay_with_same_target": pa.Column(str, required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for assays."""
