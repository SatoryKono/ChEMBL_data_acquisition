"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.common import PipelineRunMetadata, run_steps
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
    **_: object,
) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    target_columns: Sequence[str] = [
        "standard_type",
        "standard_relation",
        "standard_units",
    ]

    for column in target_columns:
        if column not in normalised.columns:
            continue

        series = normalised[column].astype("string").str.strip().str.replace(
            "\s+", " ", regex=True
        )
        if column == "standard_units" and not enforce_uppercase_units:
            normalised[column] = series
            continue
        if column in {"standard_type", "standard_relation"} and not relation_normalization:
            normalised[column] = series
            continue
        normalised[column] = series.str.upper()

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Sequence[str] | None = None,
    default_quality_flag: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    comment_column = "data_validity_comment"
    terms = [term.lower() for term in (quality_terms or ("valid",))]

    if comment_column in enriched.columns:
        comment_series = enriched[comment_column].fillna("").astype("string")

        def _quality_match(value: str) -> bool:
            lower = value.lower()
            return any(term in lower for term in terms)

        enriched["quality_flag"] = comment_series.map(_quality_match)
    else:
        enriched["quality_flag"] = default_quality_flag

    if comment_column not in enriched.columns:
        enriched[comment_column] = pd.Series(
            [pd.NA] * len(enriched), dtype="string"
        )

    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    numeric_identifier_dtype: str | type | None = "Int64",
    **_: object,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns and numeric_identifier_dtype is not None:
        prepared["activity_id"] = pd.to_numeric(
            prepared["activity_id"], errors="coerce"
        ).astype(numeric_identifier_dtype)

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
PIPELINE_STEPS = PIPELINE_CONFIG.steps


def run_activity_pipeline(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetadata]:
    """Execute the activity postprocessing steps."""

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
