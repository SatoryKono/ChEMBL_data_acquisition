"""Pandera schema describing target postprocessing outputs."""

from __future__ import annotations

from library._compat.pandera import pa

from .common import string_column

__all__ = ["TargetSchema"]


TargetSchema: pa.DataFrameSchema = pa.DataFrameSchema(
    {
        "target_chembl_id": string_column(required=True, nullable=False),
        "pref_name": string_column(required=True, nullable=False),
        "target_type": string_column(required=True, nullable=False),
        "organism": string_column(required=False, nullable=True),
        "target_class": string_column(required=False, nullable=True),
        "protein_family": string_column(required=False, nullable=True),
        "synonyms": string_column(required=False, nullable=True),
        "pipeline_version": string_column(required=False, nullable=True),
    }
)
