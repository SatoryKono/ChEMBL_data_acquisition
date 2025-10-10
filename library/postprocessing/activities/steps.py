"""Transformation steps for the activity postprocessing pipeline."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocessing.common import run_steps
from library.postprocessing.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
from library.postprocessing.common.logging import PipelineRunMetrics

from .schema import ACTIVITY_SCHEMA, validate_activities


def _normalise_string_series(
    series: pd.Series,
    *,
    uppercase: bool = False,
) -> pd.Series:
    """Return ``series`` coerced to ``StringDtype`` with whitespace collapsed."""

    normalised = (
        series.astype("string").str.strip().str.replace(r"\s+", " ", regex=True)
    )
    if uppercase:
        normalised = normalised.str.upper()
    return normalised


def normalize_activity_records(
    df: pd.DataFrame,
    *,
    relation_normalization: bool = True,
    enforce_uppercase_units: bool = True,
) -> pd.DataFrame:
    """Normalize core activity metadata columns.

    Parameters are aligned with :mod:`config.pipeline.activities` so callers can
    toggle specific coercions without editing the implementation. The defaults
    preserve the historical behaviour of trimming whitespace and uppercasing
    relation/unit fields.
    """

    normalised = df.copy(deep=True)
    normalised.columns = [col.strip().lower() for col in normalised.columns]

    if "standard_type" in normalised.columns:
        normalised["standard_type"] = _normalise_string_series(
            normalised["standard_type"], uppercase=True
        )

    if "standard_relation" in normalised.columns:
        normalised["standard_relation"] = _normalise_string_series(
            normalised["standard_relation"], uppercase=relation_normalization
        )
    if "standard_units" in normalised.columns:
        normalised["standard_units"] = _normalise_string_series(
            normalised["standard_units"], uppercase=enforce_uppercase_units
        )

    return normalised


def _prepare_quality_terms(terms: Iterable[str] | None) -> tuple[str, ...]:
    """Return lower-cased quality terms stripped of whitespace."""

    if not terms:
        return ()
    prepared: list[str] = []
    for term in terms:
        text = str(term).strip().lower()
        if text:
            prepared.append(text)
    return tuple(prepared)


def enrich_activity_quality(
    df: pd.DataFrame,
    *,
    quality_terms: Sequence[str] | None = None,
    default_quality_flag: bool = False,
) -> pd.DataFrame:
    """Add deterministic quality flags based on data validity comments."""

    enriched = df.copy(deep=True)
    prepared_terms = _prepare_quality_terms(quality_terms)
    if quality_terms is None and not prepared_terms:
        effective_terms: tuple[str, ...] = ("valid",)
    else:
        effective_terms = prepared_terms

    base_flag = pd.Series(bool(default_quality_flag), index=enriched.index)
    column = enriched.get("data_validity_comment")
    if column is None:
        enriched["quality_flag"] = base_flag.astype(bool)
        return enriched

    comment_series = column.fillna("").astype("string")
    if effective_terms:
        pattern = "|".join(re.escape(term) for term in effective_terms)
        matches = comment_series.str.contains(pattern, case=False, regex=True)
        quality_series = matches.fillna(bool(default_quality_flag)).astype(bool)
    else:
        quality_series = base_flag.astype(bool)

    enriched["quality_flag"] = quality_series
    return enriched


def finalize_activity_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    numeric_identifier_dtype: str | type[str] = "Int64",
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`ACTIVITY_SCHEMA`."""

    prepared = df.copy(deep=True)
    if "activity_id" in prepared.columns:
        numeric = pd.to_numeric(prepared["activity_id"], errors="coerce")
        try:
            prepared["activity_id"] = pd.Series(
                pd.array(numeric, dtype=numeric_identifier_dtype),
                index=prepared.index,
            )
        except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
            raise ValueError(
                "numeric_identifier_dtype must be a valid pandas nullable integer dtype"
            ) from exc

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
