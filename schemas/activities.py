"""Schema definitions for activity data.

This module provides the :data:`ActivitiesSchema` object describing
expected structure of activity dataframes.
"""

from __future__ import annotations

from pandera import Check, Column, DataFrameSchema

# Definition of the schema describing the activities table.
ActivitiesSchema: DataFrameSchema = DataFrameSchema(  # type: ignore[no-untyped-call]
    {
        "activity_id": Column(int, Check.ge(0), required=True),
        "testitem_id": Column(str, required=True),
        "target_id": Column(str, required=False),
        "standard_type": Column(
            str,
            Check.isin(["IC50", "EC50", "Ki", "Kd"]),
            required=False,
        ),
        "standard_value": Column(float, Check.ge(0), required=True),
        "pA_value": Column(float, Check.in_range(0, 14), required=False),
    }
)

"""pandera.DataFrameSchema: Validation schema for activities."""
