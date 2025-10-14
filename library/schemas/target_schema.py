"""Pandera schema for the reproducible target pipeline."""

from __future__ import annotations

import pandas as pd

from library._compat.pandera import pa

from .common import string_column

__all__ = ["TargetSchema"]


class TargetSchema:
    """Validate the canonical target export produced by the pipeline."""

    schema: pa.DataFrameSchema = pa.DataFrameSchema(
        {
            "target_chembl_id": string_column(required=True, nullable=False),
            "pref_name": string_column(required=True, nullable=False),
            "target_type": string_column(required=True, nullable=False),
            "organism": string_column(required=True, nullable=False),
            "uniprot_id": string_column(required=True, nullable=False),
            "gene_symbol": string_column(required=True, nullable=False),
            "target_class": string_column(required=False, nullable=True),
            "protein_family": string_column(required=False, nullable=True),
            "synonyms": string_column(required=False, nullable=True),
        },
        coerce=True,
        strict=True,
    )

    @classmethod
    def validate(cls, frame: pd.DataFrame) -> pd.DataFrame:
        """Return ``frame`` after Pandera validation."""

        return cls.schema.validate(frame, lazy=True)
