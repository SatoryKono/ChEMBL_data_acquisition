"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

import pandera.pandas as pa

# Definition of the schema describing the activities table.
ActivitiesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": pa.Column(object, required=True, nullable=True),
        "molecule_chembl_id": pa.Column(str, required=True, nullable=True),
        "assay_chembl_id": pa.Column(str, required=True, nullable=True),
        "activity_comment": pa.Column(str, required=False, nullable=True),
        "assay_description": pa.Column(str, required=False, nullable=True),
        "assay_variant_accession": pa.Column(str, required=False, nullable=True),
        "assay_variant_mutation": pa.Column(str, required=False, nullable=True),
        "bao_format": pa.Column(str, required=False, nullable=True),
        "bao_label": pa.Column(str, required=False, nullable=True),
        "data_validity_comment": pa.Column(str, required=False, nullable=True),
        "data_validity_description": pa.Column(str, required=False, nullable=True),
        "document_chembl_id": pa.Column(str, required=False, nullable=True),
        "pchembl_value": pa.Column(object, required=False, nullable=True),
        "potential_duplicate": pa.Column(str, required=False, nullable=True),
        "qudt_units": pa.Column(str, required=False, nullable=True),
        "record_id": pa.Column(str, required=False, nullable=True),
        "relation": pa.Column(str, required=False, nullable=True),
        "src_assay_id": pa.Column(object, required=False, nullable=True),
        "src_id": pa.Column(object, required=False, nullable=True),
        "standard_relation": pa.Column(str, required=False, nullable=True),
        "standard_units": pa.Column(str, required=False, nullable=True),
        "type": pa.Column(str, required=False, nullable=True),
        "units": pa.Column(str, required=False, nullable=True),
        "value": pa.Column(object, required=False, nullable=True),
        "standard_type": pa.Column(
            str,
            pa.Check.isin(["IC50", "Ki"]),
            required=False,
            nullable=True,
        ),
        "standard_value": pa.Column(
            object, pa.Check.ge(0), required=True, nullable=True, coerce=True
        ),
        #    "pA_value": pa.Column(float, pa.Check.in_range(-14, 14), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for activities."""
