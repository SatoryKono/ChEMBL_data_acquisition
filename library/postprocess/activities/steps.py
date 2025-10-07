"""Transformation steps for the activity postprocessing pipeline."""
from __future__ import annotations

import re
from collections.abc import Sequence

import pandas as pd

 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics

from library.pipelines.common.metadata import get_pipeline_version

from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)

from library.processing.activity.bounds import _normalize_relation as _canonical_relation


from .schema import ACTIVITY_SCHEMA, validate_activities


def normalize_activity_records(
    df: pd.DataFrame,
    *,
    relation_normalization: bool = False,
    enforce_uppercase_units: bool = True,
) -> pd.DataFrame:
    """Normalise textual activity fields and enforce consistent naming.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw activity frame to normalise.
    relation_normalization : bool, optional
        When ``True`` apply canonical mapping to ``standard_relation`` values
        (for example ``">"``/">="`` → ``">="``). Defaults to ``False`` to
        preserve legacy behaviour.
    enforce_uppercase_units : bool, optional
        When ``True`` coerce ``standard_units`` to uppercase after trimming
        whitespace. Defaults to ``True`` for backwards compatibility.
    """

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    target_columns = ["standard_type", "standard_relation", "standard_units"]
    for column in target_columns:
        if column not in normalised.columns:
            continue

        series = (
            normalised[column]
            .astype("string")
            .str.strip()
            .str.replace(r"\s+", " ", regex=True)
        )

        if column == "standard_relation" and relation_normalization:
            normalised[column] = series.map(_canonical_relation).astype("string")
        else:
            if column == "standard_units" and not enforce_uppercase_units:
                normalised[column] = series
            else:
                normalised[column] = series.str.upper()

    return normalised


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Sequence[str] | None = None,
    default_quality_flag: bool = False,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    terms = [term for term in (quality_terms or ("valid",)) if term]

    if "data_validity_comment" in enriched.columns:
        comment_series = enriched["data_validity_comment"].fillna("").astype("string")
        if terms:
            pattern = "|".join(re.escape(term) for term in terms)
            matches = comment_series.str.contains(pattern, case=False, regex=True)
        else:
            matches = pd.Series(False, index=enriched.index, dtype="boolean")
        enriched["quality_flag"] = matches.fillna(default_quality_flag)
    else:
        enriched["quality_flag"] = pd.Series(
            [default_quality_flag] * len(enriched), index=enriched.index, dtype="boolean"
        )

    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    numeric_identifier_dtype: str | type | None = "Int64",
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns and numeric_identifier_dtype is not None:
        numeric_series = pd.to_numeric(prepared["activity_id"], errors="coerce")
        try:
            prepared["activity_id"] = numeric_series.astype(numeric_identifier_dtype)
        except TypeError:
            prepared["activity_id"] = numeric_series.astype("Int64")

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
