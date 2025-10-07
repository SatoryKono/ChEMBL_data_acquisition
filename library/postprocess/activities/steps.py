"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

from library.postprocess.common import run_steps
from library.postprocess.config import load_pipeline_config

from .schema import ACTIVITY_SCHEMA, validate_activities


def normalize_activity_records(
    df: pd.DataFrame,
    *,
    strip_columns: Sequence[str] | None = None,
    uppercase_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    for column in strip_columns or ():
        if column in normalised.columns:
            normalised[column] = normalised[column].astype("string").str.strip()

    for column in uppercase_columns or ("standard_type", "standard_relation", "standard_units"):
        if column in normalised.columns:
            normalised[column] = (
                normalised[column]
                .astype("string")
                .str.strip()
                .str.replace("\s+", " ", regex=True)
                .str.upper()
            )

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    comment_column: str = "data_validity_comment",
    flag_column: str = "quality_flag",
    positive_indicators: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    indicators = {value.lower() for value in (positive_indicators or ("valid",))}
    if comment_column in enriched.columns:
        comment_series = enriched[comment_column].fillna("").astype("string")
        enriched[flag_column] = comment_series.str.lower().apply(
            lambda value: any(token in value for token in indicators)
        )
    else:
        enriched[flag_column] = False

    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    id_column: str = "activity_id",
    string_columns: Sequence[str] | None = None,
    order_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if id_column in prepared.columns:
        prepared[id_column] = pd.to_numeric(prepared[id_column], errors="coerce").astype("Int64")
    for column in string_columns or (
        "molecule_chembl_id",
        "assay_chembl_id",
        "standard_type",
        "standard_relation",
        "standard_units",
    ):
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if order_columns:
        ordered = [col for col in order_columns if col in prepared.columns]
        remaining = [col for col in prepared.columns if col not in ordered]
        prepared = prepared.loc[:, [*ordered, *remaining]]

    validated = validate_activities(prepared, context="activity_finalization")
    return validated


_PIPELINE = load_pipeline_config("activities")
PIPELINE_VERSION = _PIPELINE.pipeline_version
PIPELINE_STEPS = _PIPELINE.steps


def run_activity_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Execute the activity postprocessing steps."""

    version = pipeline_version or PIPELINE_VERSION
    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=ACTIVITY_SCHEMA,
        pipeline_version=version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "PIPELINE_VERSION",
    "finalize_activity_records",
    "normalize_activity_records",
    "run_activity_pipeline",
    "enrich_activity_quality",
]
