"""Transformation steps for assay postprocessing."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.common import PipelineRunMetadata, run_steps
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
    **_: object,
) -> pd.DataFrame:
    """Normalize string-based assay descriptors."""

    normalized = df.copy(deep=True)
    columns = ["assay_type", "assay_test_type", "assay_format"]

    for column in columns:
        if column not in normalized.columns:
            continue

        series = normalized[column].astype("string")
        if strip_whitespace:
            series = series.str.strip().str.replace("\s+", " ", regex=True)
        if uppercase_categories:
            series = series.str.upper()
        normalized[column] = series

    return normalized


def enrich_assay_flags(
    df: pd.DataFrame,
    *,
    confirmatory_terms: Sequence[str] | None = None,
    default_flag: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    type_series = enriched.get("assay_type")
    terms = [term.lower() for term in (confirmatory_terms or ("confirm",))]

    if type_series is not None:
        assay_types = type_series.astype("string").fillna("")

        def _matches(value: str) -> bool:
            lower = value.lower()
            return any(term in lower for term in terms)

        enriched["is_confirmatory"] = assay_types.map(_matches)
    else:
        enriched["is_confirmatory"] = default_flag
    return enriched


def finalize_assay_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    normalize_identifiers: bool = True,
    **_: object,
) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    if normalize_identifiers:
        for column in [
            "assay_chembl_id",
            "assay_type",
            "assay_test_type",
            "description",
        ]:
            if column in prepared.columns:
                prepared[column] = prepared[column].astype("string")

    if not enforce_schema:
        return prepared

    validated = validate_assays(prepared, context="assay_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("assays")
PIPELINE_STEPS = PIPELINE_CONFIG.steps


def run_assay_pipeline(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetadata]:
    """Run the assay postprocessing pipeline."""

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
