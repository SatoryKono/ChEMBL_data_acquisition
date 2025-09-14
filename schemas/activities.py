"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

import pandera.pandas as pa

# Definition of the schema describing the activities table.
ActivitiesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": pa.Column(str, required=True),
        "molecule_chembl_id": pa.Column(str, required=True),
        "assay_chembl_id": pa.Column(str, required=True),
        "activity_comment": pa.Column(str, required=False),
        "assay_description": pa.Column(str, required=False),
        "assay_variant_accession": pa.Column(str, required=False),
        "assay_variant_mutation": pa.Column(str, required=False),
        "bao_format": pa.Column(str, required=False),
        "bao_label": pa.Column(str, required=False),
        "data_validity_comment": pa.Column(str, required=False),
        "data_validity_description": pa.Column(str, required=False),
        "document_chembl_id": pa.Column(str, required=False),
        "pchembl_value": pa.Column(str, required=False),
        "potential_duplicate": pa.Column(str, required=False),
        "qudt_units": pa.Column(str, required=False),
        "record_id": pa.Column(str, required=False),
        "relation": pa.Column(str, required=False),
        "src_assay_id": pa.Column(str, required=False),
        "src_id": pa.Column(str, required=False),
        "standard_relation": pa.Column(str, required=False),
        "standard_units": pa.Column(str, required=False),
        "type": pa.Column(str, required=False),
        "units": pa.Column(str, required=False),
        "value": pa.Column(str, required=False),
        "standard_type": pa.Column(
            str,
            pa.Check.isin(["IC50", "Ki"]),
            required=False,
        ),
        "standard_value": pa.Column(float, pa.Check.ge(0), required=True),
        #    "pA_value": pa.Column(float, pa.Check.in_range(-14, 14), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for activities."""
