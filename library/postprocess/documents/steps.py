"""Transformation steps for document postprocessing."""
from __future__ import annotations

import pandas as pd

from library.postprocess.common import StepDefinition, run_steps

from .schema import DOCUMENT_SCHEMA, validate_documents


def normalize_document_fields(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize whitespace and casing for textual document fields."""

    normalized = df.copy(deep=True)
    normalized.columns = [col.strip().lower() for col in normalized.columns]
    for column in ["title", "journal", "doc_type"]:
        if column in normalized.columns:
            normalized[column] = normalized[column].astype("string").str.strip()
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


PIPELINE_STEPS = [
    StepDefinition("normalize_document_fields", normalize_document_fields),
    StepDefinition("enrich_document_publication_year", enrich_document_publication_year),
    StepDefinition("finalize_document_records", finalize_document_records),
]


def run_document_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Run the document postprocessing pipeline."""

    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=DOCUMENT_SCHEMA,
        pipeline_version=pipeline_version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "finalize_document_records",
    "normalize_document_fields",
    "run_document_pipeline",
    "enrich_document_publication_year",
]
