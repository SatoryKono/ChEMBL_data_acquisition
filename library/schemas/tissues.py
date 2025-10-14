"""Schema definition for tissue metadata exported by the pipeline."""

from __future__ import annotations

from typing import Any, Final, cast

from library._compat.pandera import Column, DataFrameSchema

FLEXIBLE_DTYPE: Final[Any] = cast(Any, None)

TissuesSchema: DataFrameSchema = DataFrameSchema(
    {
        "tissue_chembl_id": Column(str, required=True, nullable=True),
        "pref_name": Column(str, required=False, nullable=True),
        "uberon_id": Column(str, required=False, nullable=True),
        "efo_id": Column(str, required=False, nullable=True),
        "bto_id": Column(str, required=False, nullable=True),
        "caloha_id": Column(str, required=False, nullable=True),
        "pipeline_version": Column(str, required=False, nullable=True),
        "timestamp_utc": Column(str, required=False, nullable=True),
    }
)
"""DataFrameSchema enforcing the column layout for tissue exports."""

__all__ = ["TissuesSchema"]
