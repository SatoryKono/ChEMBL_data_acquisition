"""Transformation steps for target postprocessing."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from library.postprocess.common import run_steps
from library.postprocess.config import load_pipeline_config

from .schema import TARGET_SCHEMA, validate_targets


def normalize_target_fields(
    df: pd.DataFrame,
    *,
    id_column: str = "target_chembl_id",
    strip_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Normalize leading/trailing whitespace and ensure upper-case identifiers."""

    normalized = df.copy(deep=True)
    if id_column in normalized.columns:
        normalized[id_column] = normalized[id_column].astype("string").str.strip().str.upper()
    for column in strip_columns or ("pref_name", "organism"):
        if column in normalized.columns:
            normalized[column] = normalized[column].astype("string").str.strip()
    return normalized


def enrich_target_synonyms(
    df: pd.DataFrame,
    *,
    synonyms_column: str = "synonyms",
    separator: str = ",",
) -> pd.DataFrame:
    """Ensure synonyms column is deterministically ordered."""

    enriched = df.copy(deep=True)
    if synonyms_column in enriched.columns:
        joiner = f"{separator} " if separator and not separator.endswith(" ") else separator
        enriched[synonyms_column] = (
            enriched[synonyms_column]
            .fillna("")
            .astype("string")
            .apply(
                lambda value: joiner.join(
                    sorted(
                        part.strip()
                        for part in str(value).split(separator)
                        if part and part.strip()
                    )
                )
            )
        )
    return enriched


def finalize_target_records(
    df: pd.DataFrame,
    *,
    string_columns: Sequence[str] | None = None,
    sort_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`TARGET_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in string_columns or ("target_chembl_id", "pref_name", "target_type"):
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if sort_columns:
        existing = [col for col in sort_columns if col in prepared.columns]
        if existing:
            prepared = prepared.sort_values(existing, kind="mergesort").reset_index(drop=True)

    validated = validate_targets(prepared, context="target_finalization")
    return validated


_PIPELINE = load_pipeline_config("targets")
PIPELINE_VERSION = _PIPELINE.pipeline_version
PIPELINE_STEPS = _PIPELINE.steps


def run_target_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Run the target postprocessing pipeline."""

    version = pipeline_version or PIPELINE_VERSION
    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=TARGET_SCHEMA,
        pipeline_version=version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "PIPELINE_VERSION",
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "enrich_target_synonyms",
]
