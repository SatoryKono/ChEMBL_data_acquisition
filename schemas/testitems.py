"""Schema definitions for test item data."""

from __future__ import annotations

from pandera import Check, Column, DataFrameSchema

TestitemsSchema: DataFrameSchema = DataFrameSchema(
    {
        "salt_chembl_id": Column(str, required=True),
        "molecule_chembl_id": Column(str, required=True),
        "molecule_type": Column(
            str,
            Check.isin(["Small molecule", "Biopolymer", "Oligosaccharide", "Unknown"]),
            required=True,
        ),
        "chirality": Column(int, Check.isin([-1, 0, 1, 2]), required=False),
        "mw_freebase": Column(float, Check.in_range(0, 2000), required=True),
        "num_ro5_violations": Column(float, Check.in_range(0, 5), required=False),
        "is_radical": Column(bool, required=False),
    }
)

"""pandera.DataFrameSchema: Validation schema for test items."""
