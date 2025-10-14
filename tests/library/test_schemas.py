from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.common.schema import DataFrameSchema, validate_schema
from library.postprocessing.common.types import SchemaValidationError


@pytest.mark.unit
def test_validate_schema_coerces_types_and_preserves_rows() -> None:
    schema = DataFrameSchema(
        name="UnitTestSchema",
        required_columns=("identifier", "value"),
        dtypes={"identifier": "Int64", "value": "string"},
    )
    frame = pd.DataFrame({"identifier": ["1", "2"], "value": ["a", "b"]})

    validated = schema.validate(frame, context="unit_test")

    assert list(validated["identifier"]) == [1, 2]
    assert validated["identifier"].dtype.name == "Int64"
    assert validated["value"].dtype.name == "string"
    assert len(validated.index) == len(frame.index)


@pytest.mark.unit
def test_validate_schema_logs_failures_and_raises(capfd: pytest.CaptureFixture[str]) -> None:
    schema = DataFrameSchema(name="FailingSchema", required_columns=("identifier",))
    frame = pd.DataFrame({"value": ["missing"]})

    with pytest.raises(SchemaValidationError):
        validate_schema(frame, schema, context="unit_test")

    log_text = capfd.readouterr().err
    assert "Schema validation failed" in log_text
    assert "FailingSchema" in log_text
    assert "unit_test" in log_text
