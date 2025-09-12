"""Schema definitions for target data."""

from __future__ import annotations

import pandera.pandas as pa

TargetsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "target_chembl_id": pa.Column(str, required=True),
        # ``organism`` is not present in all datasets; treat as optional to allow
        # validation of partial records without skipping the entire step.
        "organism": pa.Column(str, required=False),
        "uniprot_id": pa.Column(str, required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for targets."""
