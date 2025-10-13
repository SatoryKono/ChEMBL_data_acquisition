"""Schema definition for cell line metadata exported by the pipeline."""

from __future__ import annotations

from typing import Any, Final, cast

from library._compat.pandera import pa

FLEXIBLE_DTYPE: Final[Any] = cast(Any, None)

CellLinesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "cell_chembl_id": pa.Column(str, required=True, nullable=True),
        "cell_name": pa.Column(str, required=False, nullable=True),
        "cell_description": pa.Column(str, required=False, nullable=True),
        "cell_id": pa.Column(int, required=False, nullable=True),
        "cell_source_organism": pa.Column(str, required=False, nullable=True),
        "cell_source_tax_id": pa.Column(int, required=False, nullable=True),
        "cell_source_tissue": pa.Column(str, required=False, nullable=True),
        "cellosaurus_id": pa.Column(str, required=False, nullable=True),
        "cl_lincs_id": pa.Column(str, required=False, nullable=True),
        "clo_id": pa.Column(str, required=False, nullable=True),
        "efo_id": pa.Column(str, required=False, nullable=True),
        "pipeline_version": pa.Column(str, required=False, nullable=True),
        "timestamp_utc": pa.Column(str, required=False, nullable=True),
    }
)
"""pandera schema enforcing the column layout for cell line exports."""

__all__ = ["CellLinesSchema"]
