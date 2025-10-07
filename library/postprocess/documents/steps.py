"""Transformation steps for document postprocessing."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from library.postprocess.common import run_steps
from library.postprocess.config import load_pipeline_config

from .schema import DOCUMENT_SCHEMA, validate_documents


def normalize_document_fields(
    df: pd.DataFrame,
    *,
    lowercase_columns: bool = True,
    strip_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Normalize whitespace and casing for textual document fields."""

    normalized = df.copy(deep=True)
    if lowercase_columns:
        normalized.columns = [col.strip().lower() for col in normalized.columns]
    else:
        normalized.columns = [col.strip() for col in normalized.columns]
    for column in strip_columns or ("title", "journal", "doc_type"):
        if column in normalized.columns:
            normalized[column] = normalized[column].astype("string").str.strip()
    return normalized


def enrich_document_publication_year(
    df: pd.DataFrame,
    *,
    source_column: str = "year",
    target_column: str = "publication_year",
) -> pd.DataFrame:
    """Create a deterministic ``publication_year`` column."""

    enriched = df.copy(deep=True)
    if source_column in enriched.columns:
        enriched[target_column] = pd.to_numeric(enriched[source_column], errors="coerce").astype(
            "Int64"
        )
    else:
        enriched[target_column] = pd.Series([pd.NA] * len(enriched), dtype="Int64")
    return enriched


def finalize_document_records(
    df: pd.DataFrame,
    *,
    string_columns: Sequence[str] | None = None,
    sort_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`DOCUMENT_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in string_columns or ("document_chembl_id", "title", "doc_type"):
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if sort_columns:
        existing = [col for col in sort_columns if col in prepared.columns]
        if existing:
            prepared = prepared.sort_values(existing, kind="mergesort").reset_index(drop=True)

    if "publication_year" in prepared.columns:
        prepared["publication_year"] = prepared["publication_year"].astype("Int64")

    validated = validate_documents(prepared, context="document_finalization")
    return validated


_PIPELINE = load_pipeline_config("documents")
PIPELINE_VERSION = _PIPELINE.pipeline_version
PIPELINE_STEPS = _PIPELINE.steps


def run_document_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Run the document postprocessing pipeline."""

    version = pipeline_version or PIPELINE_VERSION
    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=DOCUMENT_SCHEMA,
        pipeline_version=version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "PIPELINE_VERSION",
    "finalize_document_records",
    "normalize_document_fields",
    "run_document_pipeline",
    "enrich_document_publication_year",
]
