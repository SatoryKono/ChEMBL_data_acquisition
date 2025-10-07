"""Transformation steps for assay postprocessing."""
from __future__ import annotations

from typing import Iterable, Sequence

import pandas as pd

from library.postprocess.common import run_steps
from library.postprocess.config import load_pipeline_config

from .schema import ASSAY_SCHEMA, validate_assays


def normalize_assay_metadata(
    df: pd.DataFrame,
    *,
    categorical_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Normalize string-based assay descriptors."""

    normalized = df.copy(deep=True)
    for column in categorical_columns or ("assay_type", "assay_test_type", "assay_format"):
        if column in normalized.columns:
            normalized[column] = (
                normalized[column]
                .astype("string")
                .str.strip()
                .str.replace("\s+", " ", regex=True)
                .str.upper()
            )
    return normalized


def enrich_assay_flags(
    df: pd.DataFrame,
    *,
    type_column: str = "assay_type",
    flag_column: str = "is_confirmatory",
    positive_indicators: Iterable[str] | None = None,
) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    positive = {value.lower() for value in (positive_indicators or ("confirm",))}
    if type_column in enriched.columns:
        series = enriched[type_column].astype("string").str.lower()
        enriched[flag_column] = series.apply(
            lambda value: any(token in value for token in positive)
        )
    else:
        enriched[flag_column] = False
    return enriched


def finalize_assay_records(
    df: pd.DataFrame,
    *,
    string_columns: Sequence[str] | None = None,
    sort_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    for column in string_columns or (
        "assay_chembl_id",
        "assay_type",
        "assay_test_type",
        "description",
    ):
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if sort_columns:
        existing = [col for col in sort_columns if col in prepared.columns]
        if existing:
            prepared = prepared.sort_values(existing, kind="mergesort").reset_index(drop=True)

    validated = validate_assays(prepared, context="assay_finalization")
    return validated


_PIPELINE = load_pipeline_config("assays")
PIPELINE_VERSION = _PIPELINE.pipeline_version
PIPELINE_STEPS = _PIPELINE.steps


def run_assay_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> pd.DataFrame:
    """Run the assay postprocessing pipeline."""

    version = pipeline_version or PIPELINE_VERSION
    return run_steps(
        df,
        PIPELINE_STEPS,
        schema=ASSAY_SCHEMA,
        pipeline_version=version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "PIPELINE_VERSION",
    "finalize_assay_records",
    "normalize_assay_metadata",
    "run_assay_pipeline",
    "enrich_assay_flags",
]
