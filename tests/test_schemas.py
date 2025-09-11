"""Unit tests for pandera DataFrame schemas."""

from __future__ import annotations

import pandas as pd
import pytest
from pandera.errors import SchemaError

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
