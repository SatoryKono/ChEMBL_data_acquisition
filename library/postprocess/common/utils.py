"""Pipeline orchestration utilities for postprocessing transformations."""
from __future__ import annotations

from typing import Iterable, Optional

import pandas as pd

from . import logging as logging_utils
from .io import clone_dataframe, ensure_dataframe
from .schema import DataFrameSchema, validate_schema
from .types import SchemaValidationError, StepDefinition, StepError, StepIterable


def run_steps(
    df: pd.DataFrame,
    steps: StepIterable,
    *,
    schema: DataFrameSchema | None = None,
    pipeline_version: str | None = None,
    logger=None,
) -> pd.DataFrame:
    """Execute ``steps`` sequentially returning the processed DataFrame."""

    log = logger or logging_utils.get_logger()
    current = ensure_dataframe(df, copy=True)
    if pipeline_version is not None:
        current.attrs["pipeline_version"] = pipeline_version

    for step in steps:
        log.info("Starting step %s", step.name)
        before_columns = list(current.columns)
        before_shape = current.shape
        try:
            next_df = step.func(clone_dataframe(current))
        except SchemaValidationError:
            raise
        except StepError:
            raise
        except Exception as exc:  # pragma: no cover - defensive
            raise StepError(step.name, str(exc), cause=exc) from exc

        if not isinstance(next_df, pd.DataFrame):
            raise StepError(step.name, "step did not return a pandas.DataFrame")

        next_df = ensure_dataframe(next_df, copy=True)
        if pipeline_version is not None:
            next_df.attrs["pipeline_version"] = pipeline_version

        added_columns = [col for col in next_df.columns if col not in before_columns]
        removed_columns = [col for col in before_columns if col not in next_df.columns]
        log.info(
            "Completed step %s | rows=%s->%s | cols=%s->%s | +%s | -%s",
            step.name,
            before_shape[0],
            next_df.shape[0],
            before_shape[1],
            next_df.shape[1],
            ",".join(added_columns) or "-",
            ",".join(removed_columns) or "-",
        )
        current = next_df

    if schema is not None:
        log.info("Validating final schema (%s)", schema.__class__.__name__)
        current = validate_schema(current, schema, context="pipeline")
        if pipeline_version is not None:
            current.attrs["pipeline_version"] = pipeline_version

    return current


__all__ = ["run_steps"]
