"""Schema definitions for target data."""

from __future__ import annotations

import pandera.pandas as pa
from pandera import Check

TargetsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "target_chembl_id": pa.Column(str, required=True),
        "pH_dependence": pa.Column(float, required=False, checks=Check.between(0, 14)),
        "isoforms": pa.Column(float, required=False, checks=Check.ge(0)),
        # ``organism`` is not present in all datasets; treat as optional to allow
        # validation of partial records without skipping the entire step.
        "organism": pa.Column(str, required=False),
        "uniprot_id": pa.Column(str, required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for targets."""
