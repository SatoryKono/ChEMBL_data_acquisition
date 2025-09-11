"""Schema definitions for target data."""

from __future__ import annotations

import pandera.pandas as pa

TargetsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "target_chembl_id": pa.Column(str, required=True),
        "organism": pa.Column(str, required=True),
        "target_uniprot_id": pa.Column(str, required=False),
        "pH_dependence": pa.Column(float, pa.Check.in_range(0, 14), required=False),
        "isoforms": pa.Column(float, pa.Check.ge(0), required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for targets."""
