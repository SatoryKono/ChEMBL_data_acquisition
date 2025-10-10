"""Document postprocessing pipeline."""

from .schema import DOCUMENT_SCHEMA, validate_documents
from .steps import (
    PIPELINE_STEPS,
    enrich_document_publication_year,
    finalize_document_records,
    normalize_document_fields,
    run_document_pipeline,
)

__all__ = [
    "DOCUMENT_SCHEMA",
    "PIPELINE_STEPS",
    "enrich_document_publication_year",
    "finalize_document_records",
    "normalize_document_fields",
    "run_document_pipeline",
    "validate_documents",
]
