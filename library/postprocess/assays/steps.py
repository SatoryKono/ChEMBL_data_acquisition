"""Transformation steps for assay postprocessing."""

from __future__ import annotations

import re
from collections.abc import Sequence

import pandas as pd


from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
from library.postprocess.common.logging import PipelineRunMetrics

from library.pipelines.common.metadata import get_pipeline_version

from .schema import ASSAY_SCHEMA, validate_assays


_DEFAULT_CATEGORICAL_COLUMNS: tuple[str, ...] = (
    "assay_type",
    "assay_test_type",
    "assay_format",
)


def normalize_assay_metadata(
    df: pd.DataFrame,
    *,
    uppercase_categories: bool = True,
    strip_whitespace: bool = True,
    collapse_whitespace: bool = True,
    target_columns: Sequence[str] | None = None,
    **_unused: object,
) -> pd.DataFrame:
    """Normalize string-based assay descriptors.

    Parameters mirror the declarative pipeline configuration so the CLI can
    tweak behaviour without requiring code changes.  ``_unused`` captures any
    forward-compatible options while keeping the transformation pure.
    """

    normalized = df.copy(deep=True)
    columns = tuple(target_columns) if target_columns is not None else _DEFAULT_CATEGORICAL_COLUMNS
    for column in columns:
        if column not in normalized.columns:
            continue
        series = normalized[column].astype("string")
        if strip_whitespace:
            series = series.str.strip()
        if collapse_whitespace:
            series = series.str.replace(r"\s+", " ", regex=True)
        if uppercase_categories:
            series = series.str.upper()
        normalized[column] = series
    return normalized


def enrich_assay_flags(
    df: pd.DataFrame,
    *,
    confirmatory_terms: Sequence[str] | None = None,
    default_flag: bool = False,
    case_sensitive: bool = False,
    **_unused: object,
) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    type_series = enriched.get("assay_type")
    if type_series is None:
        enriched["is_confirmatory"] = pd.Series(
            bool(default_flag), index=enriched.index, dtype="bool"
        )
        return enriched

    terms = tuple(
        str(term)
        for term in (confirmatory_terms or ("confirm",))
        if str(term).strip()
    )
    if not terms:
        enriched["is_confirmatory"] = pd.Series(
            bool(default_flag), index=enriched.index, dtype="bool"
        )
        return enriched

    pattern = "|".join(re.escape(term) for term in terms)
    matches = (
        type_series.astype("string")
        .str.contains(pattern, case=case_sensitive, na=default_flag)
        .astype("bool")
    )
    enriched["is_confirmatory"] = matches
    return enriched


def finalize_assay_records(df: pd.DataFrame) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    for column in ["assay_chembl_id", "assay_type", "assay_test_type", "description"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    validated = validate_assays(prepared, context="assay_finalization")
    return validated


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
