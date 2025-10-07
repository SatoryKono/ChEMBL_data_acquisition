"""Transformation steps for document postprocessing."""
from __future__ import annotations

import unicodedata
from typing import Any

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.common import PipelineRunMetadata, run_steps
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)

from .schema import DOCUMENT_SCHEMA, validate_documents


def normalize_document_fields(
    df: pd.DataFrame,
    *,
    trim_whitespace: bool = True,
    normalise_unicode: bool = True,
    **_: object,
) -> pd.DataFrame:
    """Normalize whitespace and casing for textual document fields."""

    normalized = df.copy(deep=True)
    normalized.columns = [col.strip().lower() for col in normalized.columns]
    textual_columns = ["title", "journal", "doc_type"]

    def _normalise_value(value: pd.Series) -> pd.Series:
        series = value.astype("string")
        if trim_whitespace:
            series = series.str.strip()
        if normalise_unicode:
            series = series.map(
                lambda item: (
                    unicodedata.normalize("NFKC", item)
                    if isinstance(item, str)
                    else item
                )
            )
        return series

    for column in textual_columns:
        if column in normalized.columns:
            normalized[column] = _normalise_value(normalized[column])

    return normalized


def enrich_document_publication_year(
    df: pd.DataFrame,
    *,
    fallback_year: int | None = None,
    prefer_doi_year: bool = True,
    **_: object,
) -> pd.DataFrame:
    """Create a deterministic ``publication_year`` column."""

    enriched = df.copy(deep=True)
    if "year" in enriched.columns:
        numeric = pd.to_numeric(enriched["year"], errors="coerce").astype("Int64")
        if fallback_year is not None and not prefer_doi_year:
            numeric = numeric.fillna(fallback_year).astype("Int64")
        enriched["publication_year"] = numeric
    else:
        if fallback_year is not None and not prefer_doi_year:
            fill_value: Any = fallback_year
        else:
            fill_value = pd.NA
        enriched["publication_year"] = pd.Series(
            [fill_value] * len(enriched), dtype="Int64"
        )
    return enriched


def finalize_document_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    ensure_unique_ids: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`DOCUMENT_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in ["document_chembl_id", "title", "doc_type"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if "publication_year" in prepared.columns:
        prepared["publication_year"] = prepared["publication_year"].astype("Int64")

    if ensure_unique_ids and "document_chembl_id" in prepared.columns:
        prepared = prepared.drop_duplicates(
            subset=["document_chembl_id"], keep="first"
        ).reset_index(drop=True)

    if not enforce_schema:
        return prepared

    validated = validate_documents(prepared, context="document_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("documents")
PIPELINE_STEPS = PIPELINE_CONFIG.steps


def run_document_pipeline(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetadata]:
    """Run the document postprocessing pipeline."""

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
