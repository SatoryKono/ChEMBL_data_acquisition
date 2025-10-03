"""Document aggregation and enrichment pipelines."""

from __future__ import annotations

from . import pipeline, postprocessing, type_classifier, type_terms
from .chembl_document import get_documents

__all__ = [
  "get_documents",
  "pipeline",
  "postprocessing",
  "type_classifier",
  "type_terms",
]
