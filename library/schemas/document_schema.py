"""Pandera schema and helpers for enriched document frames.

The schema enforces the canonical ChEMBL document columns and includes
optional CrossRef/OpenAlex enrichment fields to keep metadata joins
type-safe.
"""

from __future__ import annotations

from typing import Final

import pandas as pd

from library._compat.pandera import pa

from .common import string_column

__all__ = ["DocumentSchema", "validate_document_frame"]

_BASE_COLUMNS: Final[dict[str, pa.Column]] = {
    "document_chembl_id": string_column(required=True, nullable=False),
    "doi": string_column(required=True),
    "pubmed_id": string_column(required=True),
    "title": string_column(required=True),
    "doc_type": string_column(required=True),
    "journal": string_column(required=True),
    "year": string_column(required=True),
    "doi_key": string_column(required=True),
}

_CROSSREF_COLUMNS: Final[dict[str, pa.Column]] = {
    "crossref_title": string_column(),
    "crossref_doc_type": string_column(),
    "crossref_subject": string_column(),
    "crossref_error": string_column(),
}

_OPENALEX_COLUMNS: Final[dict[str, pa.Column]] = {
    "openalex_title": string_column(),
    "openalex_doc_type": string_column(),
    "openalex_type_crossref": string_column(),
    "openalex_publication_year": string_column(),
    "openalex_error": string_column(),
}

DocumentSchema: Final[pa.DataFrameSchema] = pa.DataFrameSchema(
    {
        **_BASE_COLUMNS,
        **_CROSSREF_COLUMNS,
        **_OPENALEX_COLUMNS,
    },
    coerce=True,
    strict=False,
)
"""DataFrame schema enforcing the canonical document columns."""


def validate_document_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate ``frame`` using :data:`DocumentSchema` and return the result."""

    return DocumentSchema.validate(frame)
