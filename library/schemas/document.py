"""Pandera schema describing document postprocessing outputs."""

from __future__ import annotations

from library._compat.pandera import pa

from .common import int_column, string_column

__all__ = ["DocumentPostprocessSchema"]


DocumentPostprocessSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "document_chembl_id": string_column(required=True, nullable=False),
        "title": string_column(required=True, nullable=False),
        "doc_type": string_column(required=True, nullable=False),
        "year": string_column(required=False, nullable=True),
        "publication_year": int_column(required=False, nullable=True),
        "doi": string_column(required=False, nullable=True),
        "journal": string_column(required=False, nullable=True),
        "abstract": string_column(required=False, nullable=True),
        "pipeline_version": string_column(required=False, nullable=True),
    }
)
