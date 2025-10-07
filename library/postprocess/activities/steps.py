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

from library.processing.activity.bounds import _normalize_relation as _canonicalize_relation

from .schema import ACTIVITY_SCHEMA, validate_activities


def normalize_activity_records(
    df: pd.DataFrame,
    *,
    relation_normalization: bool = False,
    enforce_uppercase_units: bool = True,
    **_: object,
) -> pd.DataFrame:
    """Normalize string columns and enforce consistent naming."""

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    def _prepare(series: pd.Series, *, uppercase: bool) -> pd.Series:
        prepared = (
            series.astype("string")
            .str.strip()
            .str.replace("\s+", " ", regex=True)
        )
        if uppercase:
            prepared = prepared.str.upper()
        return prepared

    if "standard_type" in normalised.columns:
        normalised["standard_type"] = _prepare(
            normalised["standard_type"], uppercase=True
        )

    if "standard_units" in normalised.columns:
        normalised["standard_units"] = _prepare(
            normalised["standard_units"], uppercase=enforce_uppercase_units
        )

    if "standard_relation" in normalised.columns:
        relation_series = normalised["standard_relation"]
        if relation_normalization:
            relation_series = relation_series.map(_canonicalize_relation)
        normalised["standard_relation"] = _prepare(
            pd.Series(relation_series, index=normalised.index), uppercase=True
        )

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Sequence[str] | None = None,
    default_quality_flag: bool | None = None,
    **_: object,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    default_flag = bool(default_quality_flag) if default_quality_flag is not None else False
    terms = tuple(
        term.strip().lower()
        for term in (quality_terms or ())
        if isinstance(term, str) and term.strip()
    )

    if "data_validity_comment" in enriched.columns:
        comment_series = enriched["data_validity_comment"].fillna("").astype("string")
        if terms:
            lowered = comment_series.str.lower()
            mask = lowered.apply(lambda value: any(term in value for term in terms))
            enriched["quality_flag"] = mask.astype("boolean").fillna(default_flag)
        else:
            enriched["quality_flag"] = pd.Series(default_flag, index=enriched.index, dtype="boolean")
    else:
        enriched["quality_flag"] = pd.Series(default_flag, index=enriched.index, dtype="boolean")

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

    if enforce_schema:
        return validate_activities(prepared, context="activity_finalization")
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
