from __future__ import annotations

import pandas as pd

from library.postprocess.common import DataFrameSchema, run_steps


def _increment(df: pd.DataFrame, *, increment: int) -> pd.DataFrame:
    result = df.copy(deep=True)
    result["value"] = result["value"] + increment
    return result


def _add_flag(df: pd.DataFrame, *, flag: bool) -> pd.DataFrame:
    result = df.copy(deep=True)
    result["flag"] = flag
    return result


def test_run_steps__applies_params_and_returns_metadata() -> None:
    source = pd.DataFrame({"value": [1, 2]})

    processed, metadata = run_steps(
        source,
        [
            (_increment, {"increment": 2}),
            (_add_flag, {"flag": True}),
        ],
        pipeline_version="v1",
    )

    assert source.equals(pd.DataFrame({"value": [1, 2]}))
    assert list(processed["value"]) == [3, 4]
    assert processed["flag"].tolist() == [True, True]
    assert processed.attrs["pipeline_version"] == "v1"

    assert metadata.pipeline_version == "v1"
    assert metadata.input_rows == 2
    assert metadata.input_columns == 1
    assert metadata.output_rows == 2
    assert metadata.output_columns == 2
    assert len(metadata.steps) == 2

    first_step, second_step = metadata.steps
    assert first_step.params["increment"] == 2
    assert first_step.added_columns == ()
    assert first_step.error is None
    assert second_step.added_columns == ("flag",)
    assert second_step.columns_delta == 1
    assert second_step.error is None


def test_run_steps__validates_pre_and_post_schema() -> None:
    df = pd.DataFrame({"value": [0]})
    pre_schema = DataFrameSchema(required_columns=("value",))
    post_schema = DataFrameSchema(required_columns=("value", "flag"))

    processed, metadata = run_steps(
        df,
        [
            (_increment, {"increment": 1}),
            (_add_flag, {"flag": False}),
        ],
        pipeline_version="pipeline-1",
        pre_schema=pre_schema,
        post_schema=post_schema,
    )

    assert list(processed.columns) == ["value", "flag"]
    assert processed.attrs["pipeline_version"] == "pipeline-1"
    assert metadata.pre_validation is not None
    assert metadata.pre_validation.error is None
    assert metadata.post_validation is not None
    assert metadata.post_validation.error is None
