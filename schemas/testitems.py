"""Schema definitions for test item data."""

from __future__ import annotations

import pandera.pandas as pa

TestitemsSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "salt_chembl_id": pa.Column(str, required=True),
        "salt_chembl_id": pa.Column(str, required=False),
        "molecule_chembl_id": pa.Column(str, required=True),
        "molecule_type": pa.Column(
            str,
            pa.Check.isin(
                ["Small molecule", "Biopolymer", "Oligosaccharide", "Unknown"]
            ),
            required=True,
        ),
        "chirality": pa.Column(int, pa.Check.isin([-1, 0, 1, 2]), required=False),
        "mw_freebase": pa.Column(float, pa.Check.in_range(0, 2000), required=False),
        "num_ro5_violations": pa.Column(float, pa.Check.in_range(0, 5), required=False),
        "is_radical": pa.Column(bool, required=False),
    }
)

"""pa.DataFrameSchema: Validation schema for test items."""
