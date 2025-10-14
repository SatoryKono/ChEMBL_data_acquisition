"""Schema definition for cell line metadata exported by the pipeline."""

from __future__ import annotations

from typing import Any, Final, cast

from library._compat.pandera import Column, DataFrameSchema

FLEXIBLE_DTYPE: Final[Any] = cast(Any, None)

CellLinesSchema: DataFrameSchema = DataFrameSchema(
    {
        "cell_chembl_id": Column(str, required=True, nullable=True),
        "cell_name": Column(str, required=False, nullable=True),
        "cell_description": Column(str, required=False, nullable=True),
        "cell_id": Column(int, required=False, nullable=True),
        "cell_source_organism": Column(str, required=False, nullable=True),
        "cell_source_tax_id": Column(int, required=False, nullable=True),
        "cell_source_tissue": Column(str, required=False, nullable=True),
        "cellosaurus_id": Column(str, required=False, nullable=True),
        "cl_lincs_id": Column(str, required=False, nullable=True),
        "clo_id": Column(str, required=False, nullable=True),
        "efo_id": Column(str, required=False, nullable=True),
        "pipeline_version": Column(str, required=False, nullable=True),
        "timestamp_utc": Column(str, required=False, nullable=True),
    }
)
"""DataFrameSchema enforcing the column layout for cell line exports."""

__all__ = ["CellLinesSchema"]
