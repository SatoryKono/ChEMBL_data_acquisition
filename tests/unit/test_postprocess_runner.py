from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.common import (
    DataFrameSchema,
    RunnerResult,
    StepDefinition,
    run_steps,
)


@pytest.mark.unit
def test_run_steps__applies_steps_and_records_metadata() -> None:
    frame = pd.DataFrame({"value": [1, 2]})

    schema = DataFrameSchema(required_columns=("value",), optional_columns=("total",))

    def add_constant(df: pd.DataFrame, increment: int) -> pd.DataFrame:
        result = df.copy()
        result["value"] = result["value"] + increment
        return result

    def add_total(df: pd.DataFrame, source: str) -> pd.DataFrame:
        result = df.copy()
        result["total"] = result[source].cumsum()
        return result

    result = run_steps(
        frame,
        [
            StepDefinition("add_constant", add_constant, parameters={"increment": 2}),
            ("add_total", add_total, {"source": "value"}),
        ],
        post_schema=schema,
        pipeline_version="1.2.3",
    )

    assert isinstance(result, RunnerResult)
    assert result.frame.attrs["pipeline_version"] == "1.2.3"
    assert result.metadata.pipeline_version == "1.2.3"
    assert result.metadata.pre_schema is None
    assert result.metadata.post_schema is not None
    assert result.metadata.post_schema.status == "success"

    assert result.frame["value"].tolist() == [3, 4]
    assert result.frame["total"].tolist() == [3, 7]

    step_one, step_two = result.metadata.steps
    assert step_one.name == "add_constant"
    assert step_one.parameters["increment"] == 2
    assert step_one.row_delta == 0
    assert step_one.column_delta == 0
    assert step_one.error is None

    assert step_two.name == "add_total"
    assert step_two.parameters["source"] == "value"
    assert step_two.column_delta == 1
    assert step_two.added_columns == ("total",)
    assert step_two.error is None


@pytest.mark.unit
def test_run_steps__validates_input_schema() -> None:
    frame = pd.DataFrame({"value": [1]})

    schema = DataFrameSchema(required_columns=("value",))

    def identity(df: pd.DataFrame) -> pd.DataFrame:
        return df

    result = run_steps(
        frame,
        [identity],
        pre_schema=schema,
        post_schema=schema,
        pipeline_version="0.0.1",
    )

    assert result.metadata.pre_schema is not None
    assert result.metadata.pre_schema.status == "success"
    assert result.metadata.post_schema is not None
    assert result.metadata.post_schema.status == "success"
    assert result.frame.attrs["pipeline_version"] == "0.0.1"
