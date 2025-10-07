"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

from typing import Iterable

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
    relation_normalization: bool = False,
    enforce_uppercase_units: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    def _clean_column(series: pd.Series, *, uppercase: bool) -> pd.Series:
        cleaned = series.astype("string").str.strip().str.replace("\s+", " ", regex=True)
        if uppercase:
            cleaned = cleaned.str.upper()
        return cleaned

    if "standard_type" in normalised.columns:
        normalised["standard_type"] = _clean_column(
            normalised["standard_type"], uppercase=True
        )

    if "standard_relation" in normalised.columns:
        relation_series = _clean_column(
            normalised["standard_relation"], uppercase=True
        )
        if relation_normalization:
            mapping = {
                "=": "=",
                "==": "=",
                "EQ": "=",
                "EQUAL": "=",
                "EQUALS": "=",
                "=~": "=",
                "<": "<",
                "LT": "<",
                "LESS THAN": "<",
                "<=": "<=",
                "LE": "<=",
                "LESS THAN OR EQUAL": "<=",
                "≤": "<=",
                ">": ">",
                "GT": ">",
                "GREATER THAN": ">",
                ">=": ">=",
                "GE": ">=",
                "GREATER THAN OR EQUAL": ">=",
                "≥": ">=",
                "~": "~",
                "≈": "~",
                "APPROX": "~",
                "ABOUT": "~",
            }
            relation_series = relation_series.map(
                lambda value: mapping.get(value, value)
            )
        normalised["standard_relation"] = relation_series

    if "standard_units" in normalised.columns:
        normalised["standard_units"] = _clean_column(
            normalised["standard_units"], uppercase=enforce_uppercase_units
        )

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Iterable[str] | None = None,
    default_quality_flag: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    if "data_validity_comment" in enriched.columns:
        comment_series = enriched["data_validity_comment"].fillna("").astype("string")
        terms = {
            term.strip().lower()
            for term in (quality_terms or [])
            if isinstance(term, str) and term.strip()
        }
        if not terms:
            terms = {"valid"}
        enriched["quality_flag"] = comment_series.str.lower().apply(
            lambda value: any(keyword in value for keyword in terms)
        )
    else:
        enriched["quality_flag"] = default_quality_flag

    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    numeric_identifier_dtype: str = "Int64",
    **_: object,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns:
        coerced = pd.to_numeric(prepared["activity_id"], errors="coerce")
        try:
            prepared["activity_id"] = coerced.astype(numeric_identifier_dtype)
        except TypeError:
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

    if enforce_schema:
        prepared = validate_activities(prepared, context="activity_finalization")
    return prepared


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
