"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.common import run_steps
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)

from .schema import ACTIVITY_SCHEMA, validate_activities


def normalize_activity_records(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    for column in ["standard_type", "standard_relation", "standard_units"]:
        if column in normalised.columns:
            normalised[column] = (
                normalised[column]
                .astype("string")
                .str.strip()
                .str.replace("\s+", " ", regex=True)
                .str.upper()
            )

    return normalised


def enrich_activity_quality(df: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    if "data_validity_comment" in enriched.columns:
        comment_series = enriched["data_validity_comment"].fillna("").astype("string")
        enriched["quality_flag"] = comment_series.str.contains("valid", case=False)
    else:
        enriched["quality_flag"] = False

    return enriched


def finalize_activity_records(df: pd.DataFrame) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns:
        prepared["activity_id"] = pd.to_numeric(prepared["activity_id"], errors="coerce").astype(
            "Int64"
        )
    for column in ["molecule_chembl_id", "assay_chembl_id", "standard_type", "standard_relation", "standard_units"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    validated = validate_activities(prepared, context="activity_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("activities")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_activity_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Execute the activity postprocessing steps."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=ACTIVITY_SCHEMA,
        pipeline_version=resolved_version,
        logger=logger,
    )


def _resolve_pipeline_version(override: str | None) -> str:
    """Resolve the effective pipeline version for activity postprocessing."""

    candidate = normalize_pipeline_version(override)
    if candidate is not None:
        return candidate

    config_candidate = normalize_pipeline_version(PIPELINE_CONFIG.pipeline_version)
    if config_candidate is not None:
        return config_candidate

    return get_pipeline_version()


__all__ = [
    "PIPELINE_CONFIG",
    "PIPELINE_STEPS",
    "finalize_activity_records",
    "normalize_activity_records",
    "run_activity_pipeline",
    "enrich_activity_quality",
]
