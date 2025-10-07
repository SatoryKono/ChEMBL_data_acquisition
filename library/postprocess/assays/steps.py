"""Transformation steps for assay postprocessing."""
from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import ASSAY_SCHEMA, validate_assays


def normalize_assay_metadata(
    df: pd.DataFrame,
    *,
    uppercase_categories: bool = True,
    strip_whitespace: bool = True,
    collapse_internal_whitespace: bool = True,
) -> pd.DataFrame:
    """Normalize string-based assay descriptors."""

    normalized = df.copy(deep=True)
    for column in ["assay_type", "assay_test_type", "assay_format"]:
        if column not in normalized.columns:
            continue
        series = normalized[column].astype("string")
        if strip_whitespace:
            series = series.str.strip()
        if collapse_internal_whitespace:
            series = series.str.replace(r"\s+", " ", regex=True)
        if uppercase_categories:
            series = series.str.upper()
        normalized[column] = series
    return normalized


def enrich_assay_flags(
    df: pd.DataFrame,
    *,
    confirmatory_terms: Iterable[str] | None = None,
    default_flag: bool = False,
) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    type_series = enriched.get("assay_type")
    if type_series is None:
        enriched["is_confirmatory"] = bool(default_flag)
        return enriched

    terms = tuple(term for term in (confirmatory_terms or ()) if term)
    enriched["is_confirmatory"] = bool(default_flag)
    if not terms:
        return enriched

    pattern = "|".join(re.escape(term) for term in terms)
    matches = type_series.astype("string").str.contains(
        pattern,
        case=False,
        na=False,
    )
    enriched.loc[matches, "is_confirmatory"] = True
    return enriched


def finalize_assay_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    normalize_identifiers: bool = True,
) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    if normalize_identifiers:
        for column in ["assay_chembl_id", "assay_type", "assay_test_type", "description"]:
            if column in prepared.columns:
                prepared[column] = prepared[column].astype("string")

    if enforce_schema:
        return validate_assays(prepared, context="assay_finalization")
    return prepared


PIPELINE_CONFIG = load_pipeline_config("assays")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_assay_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the assay postprocessing pipeline and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=ASSAY_SCHEMA,
        pipeline_version=resolved_version,
        logger=logger,
    )


def _resolve_pipeline_version(override: str | None) -> str:
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
    "finalize_assay_records",
    "normalize_assay_metadata",
    "run_assay_pipeline",
    "enrich_assay_flags",
]
