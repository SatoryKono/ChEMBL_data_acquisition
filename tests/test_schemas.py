"""Unit tests for pandera DataFrame schemas."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import strategies as st
from hypothesis.extra.pandas import column, data_frames, range_indexes
from pandera.errors import SchemaError, SchemaErrors

from library.normalization import normalize_activities
from schemas import (
    ActivitiesSchema,
    AssaysSchema,
    DocumentsSchema,
    TargetsSchema,
    TestitemsSchema,
)


def test_activities_schema_validation() -> None:
    """Ensure :data:`ActivitiesSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "activity_id": [1],
            "testitem_id": ["CHEMBL1"],
            "target_id": ["CHEMBL2"],
            "standard_type": ["IC50"],
            "standard_value": [10.0],
            "pA_value": [5.0],
        }
    )
    ActivitiesSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "standard_type"] = "Unknown"
    with pytest.raises(SchemaError):
        ActivitiesSchema.validate(invalid)


def test_assays_schema_validation() -> None:
    """Ensure :data:`AssaysSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1"],
            "document_chembl_id": ["CHEMBL2"],
            "target_chembl_id": ["CHEMBL3"],
            "year": [2020],
            "month": [6],
        }
    )
    AssaysSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "month"] = 13
    with pytest.raises(SchemaError):
        AssaysSchema.validate(invalid)


def test_documents_schema_validation() -> None:
    """Ensure :data:`DocumentsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1"],
            "doi": ["10.1000/xyz123"],
            "title": ["Example"],
            "year": [2020],
            "month": [1],
            "day": [15],
            "citation": [0],
        }
    )
    DocumentsSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "citation"] = -1
    with pytest.raises(SchemaError):
        DocumentsSchema.validate(invalid)


def test_targets_schema_validation() -> None:
    """Ensure :data:`TargetsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "organism": ["Homo sapiens"],
            "target_uniprot_id": ["P12345"],
            "pH_dependence": [7.0],
            "isoforms": [2.0],
        }
    )
    TargetsSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "pH_dependence"] = 20.0
    with pytest.raises(SchemaError):
        TargetsSchema.validate(invalid)


def test_testitems_schema_validation() -> None:
    """Ensure :data:`TestitemsSchema` validates expected data."""
    valid = pd.DataFrame(
        {
            "salt_chembl_id": ["CHEMBL1"],
            "molecule_chembl_id": ["CHEMBL1"],
            "molecule_type": ["Small molecule"],
            "chirality": [1],
            "mw_freebase": [100.0],
            "num_ro5_violations": [0.0],
            "is_radical": [False],
        }
    )
    TestitemsSchema.validate(valid)

    invalid = valid.copy()
    invalid.loc[0, "molecule_type"] = "Peptide"
    with pytest.raises(SchemaError):
        TestitemsSchema.validate(invalid)


@given(
    data_frames(
        columns=[
            column(
                "activity_id",
                dtype=int,
                elements=st.integers(min_value=0, max_value=10),
            ),
            column("testitem_id", elements=st.text(min_size=1)),
            column(
                "standard_value",
                elements=st.floats(min_value=0, allow_nan=False),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_activities_schema_hypothesis_valid(df: pd.DataFrame) -> None:
    """Random valid frames pass ``ActivitiesSchema``."""

    ActivitiesSchema.validate(df)


@given(
    data_frames(
        columns=[
            column(
                "activity_id",
                dtype=int,
                elements=st.integers(min_value=0, max_value=10),
            ),
            column("testitem_id", elements=st.text(min_size=1)),
            column(
                "standard_value",
                elements=st.floats(max_value=-1.0, allow_nan=False),
            ),
        ],
        index=range_indexes(min_size=1, max_size=5),
    )
)
def test_activities_schema_hypothesis_invalid(df: pd.DataFrame) -> None:
    """Negative ``standard_value`` fails validation."""

    with pytest.raises(SchemaError):
        ActivitiesSchema.validate(df)


def test_activities_from_files() -> None:
    """Validate activities data from CSV/JSON files and capture failures."""
    data_dir = Path(__file__).parent / "data"

    # Positive case: CSV data passes validation after normalisation
    valid = pd.read_csv(data_dir / "activities_valid.csv")
    normalized = normalize_activities(valid)
    assert normalized.loc[0, "relation"] == "<="
    assert normalized.loc[0, "units"] == "5 uM"
    ActivitiesSchema.validate(normalized)

    # Negative case: JSON data with invalid standard_value
    invalid = pd.read_json(data_dir / "activities_invalid.json")
    invalid_norm = normalize_activities(invalid)
    with pytest.raises(SchemaErrors) as exc_info:
        ActivitiesSchema.validate(invalid_norm, lazy=True)

    failure_cases = exc_info.value.failure_cases
    assert (failure_cases["column"] == "standard_value").any()
    assert (failure_cases["index"] == 0).any()
