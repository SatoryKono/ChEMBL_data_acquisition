"""Validation helpers for input data."""
from __future__ import annotations

from typing import Iterable

import pandas as pd

REQUIRED_COLUMNS = {
    "record_id",
    "doi",
    "pmid",
    "pmcid",
    "openalex_id",
    "pubmed_publication_types",
    "crossref_type",
    "openalex_type",
    "scholar_type",
    "pubmed_mesh_descriptors",
    "pubmed_mesh_qualifiers",
    "openalex_mesh_descriptors",
    "openalex_mesh_qualifiers",
}


def ensure_columns(df: pd.DataFrame) -> None:
    """Raise ``ValueError`` if required columns are missing."""

    missing = REQUIRED_COLUMNS.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
