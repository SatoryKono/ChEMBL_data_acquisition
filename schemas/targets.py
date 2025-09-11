"""Schema definitions for target data."""

from __future__ import annotations

from pandera import Check, Column, DataFrameSchema


TargetsSchema: DataFrameSchema = DataFrameSchema(
    {
        "target_chembl_id": Column(str, required=True),
        "organism": Column(str, required=True),
        "target_uniprot_id": Column(str, required=False),
        "pH_dependence": Column(float, Check.in_range(0, 14), required=False),
        "isoforms": Column(float, Check.ge(0), required=False),
    }
)

"""pandera.DataFrameSchema: Validation schema for targets."""
