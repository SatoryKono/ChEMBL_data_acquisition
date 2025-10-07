"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import ACTIVITY_SCHEMA, validate_activities


def normalize_activity_records(
    df: pd.DataFrame,
    *,
    relation_normalization: bool = True,
    enforce_uppercase_units: bool = True,
) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    def _clean_column(column: str, *, uppercase: bool) -> None:
        if column not in normalised.columns:
            return
        series = (
            normalised[column]
            .astype("string")
            .str.strip()
            .str.replace(r"\s+", " ", regex=True)
        )
        if uppercase:
            series = series.str.upper()
        normalised[column] = series

    _clean_column("standard_type", uppercase=True)
    _clean_column("standard_relation", uppercase=relation_normalization)
    _clean_column("standard_units", uppercase=enforce_uppercase_units)

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Sequence[str] | None = None,
    default_quality_flag: bool = False,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    if "data_validity_comment" in enriched.columns:
        comment_series = enriched["data_validity_comment"].fillna("").astype("string")
        lowered = comment_series.str.lower()
        terms = [str(term).strip().lower() for term in quality_terms or () if str(term).strip()]
        if not terms:
            terms = ["valid"]
        matches = pd.Series(False, index=lowered.index, dtype="boolean")
        for term in terms:
            matches = matches | lowered.str.contains(term, na=False)
        enriched["quality_flag"] = matches.fillna(default_quality_flag)
    else:
        enriched["quality_flag"] = bool(default_quality_flag)

    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    numeric_identifier_dtype: str = "Int64",
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns:
        coerced = pd.to_numeric(prepared["activity_id"], errors="coerce")
        try:
            prepared["activity_id"] = coerced.astype(numeric_identifier_dtype)
        except (TypeError, ValueError):  # pragma: no cover - defensive
            prepared["activity_id"] = coerced.astype("Int64")
    for column in [
        "molecule_chembl_id",
        "assay_chembl_id",
        "standard_type",
        "standard_relation",
        "standard_units",
    ]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if not enforce_schema:
        return prepared

    validated = validate_activities(prepared, context="activity_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("activities")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_activity_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Execute the activity postprocessing steps and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=ACTIVITY_SCHEMA,
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
