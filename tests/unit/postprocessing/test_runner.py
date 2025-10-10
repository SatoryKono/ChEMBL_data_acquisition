from __future__ import annotations

import logging

import pandas as pd
import pytest

from library.postprocess.common import StepDefinition, run_steps, runner
from library.postprocess.common.schema import DataFrameSchema
from library.postprocess.common.types import SchemaValidationError


def _with_constant(
    df: pd.DataFrame, *, column: str, value: int
) -> pd.DataFrame:  # pragma: no cover - exercised via tests
    result = df.copy(deep=True)
    result[column] = value
    return result


def _add_columns(
    df: pd.DataFrame, *, target: str, left: str, right: str
) -> pd.DataFrame:  # pragma: no cover - exercised via tests
    result = df.copy(deep=True)
    result[target] = result[left] + result[right]
    return result


def test_run_steps__applies_parameters_and_records_metadata() -> None:
    source = pd.DataFrame({"seed": [1, 2]}, dtype="int64")

    steps = (
        StepDefinition(
            name="with_constant",
            func=_with_constant,
            params={"column": "value", "value": 3},
        ),
        StepDefinition(
            name="combine_columns",
            func=_add_columns,
            params={"target": "total", "left": "seed", "right": "value"},
        ),
    )

    schema = DataFrameSchema(
        required_columns=("seed", "value", "total"),
        dtypes={"seed": "int64", "value": "int64", "total": "int64"},
    )

    frame, metadata = run_steps(
        source,
        steps,
        pipeline_version="2024.1",
        post_schema=schema,
    )

    assert list(frame.columns) == ["seed", "value", "total"]
    assert frame["value"].tolist() == [3, 3]
    assert frame["total"].tolist() == [4, 5]

    assert "value" not in source.columns  # original frame remains untouched
    assert "pipeline_version" not in source.attrs

    assert frame.attrs["pipeline_version"] == "2024.1"
    assert metadata.pipeline_version == "2024.1"

    assert len(metadata.steps) == 2
    assert metadata.steps[0].parameters == {"column": "value", "value": 3}
    assert metadata.steps[0].error is None
    assert metadata.steps[1].parameters["target"] == "total"


def test_run_steps__validates_input_schema_before_steps() -> None:
    frame = pd.DataFrame({"unexpected": [1]}, dtype="int64")
    pre_schema = DataFrameSchema(required_columns=("seed",))

    with pytest.raises(SchemaValidationError):
        run_steps(
            frame,
            (),
            pipeline_version="2024.1",
            pre_schema=pre_schema,
        )


def test_run_steps__accepts_tuple_step_specification() -> None:
    frame = pd.DataFrame({"seed": [0, 1]}, dtype="int64")

    final, metadata = run_steps(
        frame,
        ((_with_constant, {"column": "value", "value": 10}),),
        pipeline_version="2024.2",
        post_schema=DataFrameSchema(required_columns=("seed", "value")),
    )

    assert final["value"].tolist() == [10, 10]
    assert metadata.steps[0].name == "_with_constant"
    assert metadata.steps[0].parameters == {"column": "value", "value": 10}


def test_run_steps__ignores_unsupported_parameters(
    caplog: pytest.LogCaptureFixture,
) -> None:
    frame = pd.DataFrame({"seed": [5]}, dtype="int64")

    steps = (
        StepDefinition(
            name="with_constant",
            func=_with_constant,
            params={"column": "value", "value": 7, "unexpected": "ignored"},
        ),
    )

    test_logger = logging.getLogger("tests.postprocess.runner")
    test_logger.handlers.clear()
    test_logger.propagate = True

    with caplog.at_level("WARNING", logger="tests.postprocess.runner"):
        result, metadata = run_steps(frame, steps, logger=test_logger)

    assert result["value"].tolist() == [7]
    assert any(
        "ignoring unsupported parameters" in record.message for record in caplog.records
    )
    assert metadata.steps[0].parameters["unexpected"] == "ignored"


def test_run_steps__recovers_from_runtime_parameter_mismatch(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    frame = pd.DataFrame({"seed": [5]}, dtype="int64")

    steps = (
        StepDefinition(
            name="with_constant",
            func=_with_constant,
            params={"column": "value", "value": 11, "unexpected": "ignored"},
        ),
    )

    original_signature = runner.inspect.signature

    def _raise_signature(func):
        if func is _with_constant:
            raise ValueError("no signature available")
        return original_signature(func)

    monkeypatch.setattr(runner.inspect, "signature", _raise_signature)

    test_logger = logging.getLogger("tests.postprocess.runner")
    test_logger.handlers.clear()
    test_logger.propagate = True

    with caplog.at_level("WARNING", logger="tests.postprocess.runner"):
        result, metadata = run_steps(frame, steps, logger=test_logger)

    assert result["value"].tolist() == [11]
    assert any(
        "retrying without unsupported parameter" in record.message
        for record in caplog.records
    )
    assert metadata.steps[0].parameters["unexpected"] == "ignored"


def test_run_steps__handles_unparsed_typeerror_messages(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    frame = pd.DataFrame({"seed": [5]}, dtype="int64")

    steps = (
        StepDefinition(
            name="with_constant",
            func=_with_constant,
            params={"column": "value", "value": 13, "unexpected": "ignored"},
        ),
    )

    def _no_keyword(_error):
        return None

    def _pass_through_prepare(func, params, *, logger, step_name):
        return dict(params)

    monkeypatch.setattr(runner, "_prepare_step_arguments", _pass_through_prepare)
    monkeypatch.setattr(runner, "_extract_unexpected_keyword", _no_keyword)

    test_logger = logging.getLogger("tests.postprocess.runner")
    test_logger.handlers.clear()
    test_logger.propagate = True

    with caplog.at_level("WARNING", logger="tests.postprocess.runner"):
        result, metadata = run_steps(frame, steps, logger=test_logger)

    assert result["value"].tolist() == [13]
    assert any(
        "retrying without unsupported parameters" in record.message
        for record in caplog.records
    )
    assert metadata.steps[0].parameters["unexpected"] == "ignored"
