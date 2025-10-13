"""Schema definition for tissue metadata exported by the pipeline."""

from __future__ import annotations

from typing import Any, Final, cast

from library._compat.pandera import pa

FLEXIBLE_DTYPE: Final[Any] = cast(Any, None)

TissuesSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "tissue_chembl_id": pa.Column(str, required=True, nullable=True),
        "pref_name": pa.Column(str, required=False, nullable=True),
        "uberon_id": pa.Column(str, required=False, nullable=True),
        "efo_id": pa.Column(str, required=False, nullable=True),
        "bto_id": pa.Column(str, required=False, nullable=True),
        "caloha_id": pa.Column(str, required=False, nullable=True),
        "pipeline_version": pa.Column(str, required=False, nullable=True),
        "timestamp_utc": pa.Column(str, required=False, nullable=True),
    }
)
"""pandera schema enforcing the column layout for tissue exports."""

__all__ = ["TissuesSchema"]
