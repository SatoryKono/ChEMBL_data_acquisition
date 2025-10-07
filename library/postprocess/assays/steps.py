"""Transformation steps for assay postprocessing."""
from __future__ import annotations

import pandas as pd

 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import ASSAY_SCHEMA, validate_assays


def normalize_assay_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize string-based assay descriptors."""

    normalized = df.copy(deep=True)
    for column in ["assay_type", "assay_test_type", "assay_format"]:
        if column in normalized.columns:
            normalized[column] = (
                normalized[column]
                .astype("string")
                .str.strip()
                .str.replace("\s+", " ", regex=True)
                .str.upper()
            )
    return normalized


def enrich_assay_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    type_series = enriched.get("assay_type")
    if type_series is not None:
        enriched["is_confirmatory"] = type_series.astype("string").str.contains(
            "CONFIRM", case=False, na=False
        )
    else:
        enriched["is_confirmatory"] = False
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
