"""Transformation steps for document postprocessing."""
from __future__ import annotations

import unicodedata

import pandas as pd
 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import DOCUMENT_SCHEMA, validate_documents


def normalize_document_fields(
    df: pd.DataFrame,
    *,
    trim_whitespace: bool = True,
    normalise_unicode: bool = False,
) -> pd.DataFrame:
    """Normalize whitespace and casing for textual document fields."""

    normalized = df.copy(deep=True)
    normalized.columns = [col.strip().lower() for col in normalized.columns]

    textual_columns = ["title", "journal", "doc_type"]
    for column in textual_columns:
        if column not in normalized.columns:
            continue
        series = normalized[column].astype("string")
        if trim_whitespace:
            series = series.str.strip().str.replace(r"\s+", " ", regex=True)
        if normalise_unicode:
            series = series.map(
                lambda value: unicodedata.normalize("NFKC", value)
                if not pd.isna(value)
                else value
            )
        normalized[column] = series.replace({"": pd.NA})

    return normalized


def enrich_document_publication_year(df: pd.DataFrame) -> pd.DataFrame:
    """Create a deterministic ``publication_year`` column."""

    enriched = df.copy(deep=True)
    if "year" in enriched.columns:
        enriched["publication_year"] = pd.to_numeric(
            enriched["year"], errors="coerce"
        ).astype("Int64")
    else:
        enriched["publication_year"] = pd.Series([pd.NA] * len(enriched), dtype="Int64")
    return enriched


def finalize_document_records(df: pd.DataFrame) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`DOCUMENT_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in ["document_chembl_id", "title", "doc_type"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if "publication_year" in prepared.columns:
        prepared["publication_year"] = prepared["publication_year"].astype("Int64")

    validated = validate_documents(prepared, context="document_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("documents")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_document_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the document postprocessing pipeline and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=DOCUMENT_SCHEMA,
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
    "finalize_document_records",
    "normalize_document_fields",
    "run_document_pipeline",
    "enrich_document_publication_year",
]
