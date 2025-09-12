"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

import pandera.pandas as pa

# Definition of the schema describing the activities table.
ActivitiesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "activity_id": pa.Column(int, pa.Check.ge(0), required=True),
        "molecule_chembl_id": pa.Column(str, required=True),
        "target_id": pa.Column(str, required=False),
        "standard_type": pa.Column(
            str,
            pa.Check.isin(["IC50", "EC50", "Ki", "Kd"]),
            required=False,
        ),
        "standard_value": pa.Column(float, pa.Check.ge(0), required=True),
        "pA_value": pa.Column(float, pa.Check.in_range(0, 14), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for activities."""
