"""Schema definitions for assay post-processing."""

from __future__ import annotations

from library._compat.pandera import Column, DataFrameSchema

AssayPostprocessSchema: DataFrameSchema = DataFrameSchema(
    {
        "document_chembl_id": Column(str, required=True),
        "target_chembl_id": Column(str, required=True),
    }
)
"""DataFrameSchema: Required columns for assay post-processing."""
