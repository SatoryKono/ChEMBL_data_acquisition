"""Unit tests for postprocessing schema validation helpers."""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from library.postprocessing.activities.schema import ACTIVITY_SCHEMA
from library.postprocessing.assays.schema import ASSAY_SCHEMA
from library.postprocessing.common.schema import validate_schema
from library.postprocessing.common.types import SchemaValidationError


def test_validate_schema_success_logs(caplog: pytest.LogCaptureFixture) -> None:
    """Successful validation should emit an info log entry and return a frame."""

    frame = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL1"],
            "assay_type": ["BINDING"],
            "assay_test_type": ["PRIMARY"],
            "description": ["example"],
            "confidence_score": [0.9],
        }
    )

    caplog.set_level(logging.INFO)
    validated = validate_schema(frame, ASSAY_SCHEMA, context="assay_unit_test")

    assert list(validated.columns[:4]) == [
        "assay_chembl_id",
        "assay_type",
        "assay_test_type",
        "description",
    ]
    assert "Schema validation succeeded" in caplog.text


def test_validate_schema_failure_raises(caplog: pytest.LogCaptureFixture) -> None:
    """Pandera validation errors must propagate as SchemaValidationError."""

    frame = pd.DataFrame(
        {
            "activity_id": [1],
            "molecule_chembl_id": ["CHEMBL2"],
            "assay_chembl_id": ["CHEMBL3"],
            "standard_type": ["IC50"],
            "standard_relation": ["=",],
            "standard_value": ["invalid"],
            "standard_units": ["nM"],
        }
    )

    caplog.set_level(logging.ERROR)
    with pytest.raises(SchemaValidationError) as excinfo:
        validate_schema(frame, ACTIVITY_SCHEMA, context="activity_unit_test")

    assert "failure" in str(excinfo.value).lower()
    assert "Schema validation failed" in caplog.text
