"""Schema definitions for assay post-processing."""

from __future__ import annotations

from library._compat.pandera import pa

AssayPostprocessSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "document_chembl_id": pa.Column(str, required=True),
        "target_chembl_id": pa.Column(str, required=True),
    }
)
"""pa.DataFrameSchema: Required columns for assay post-processing."""
