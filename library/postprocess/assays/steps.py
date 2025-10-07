"""Transformation steps for assay postprocessing."""
from __future__ import annotations

import pandas as pd

from library.postprocess.common import RunnerResult, StepDefinition, run_steps

from .schema import ASSAY_SCHEMA, validate_assays


def normalize_assay_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize string-based assay descriptors."""

    normalized = df.copy(deep=True)
    for column in ["assay_type", "assay_test_type", "assay_format"]:
        if column in normalized.columns:
            normalized[column] = (
                normalized[column]
                .astype("string")
                .str.strip()
                .str.replace("\s+", " ", regex=True)
                .str.upper()
            )
    return normalized


def enrich_assay_flags(df: pd.DataFrame) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    type_series = enriched.get("assay_type")
    if type_series is not None:
        enriched["is_confirmatory"] = type_series.astype("string").str.contains(
            "CONFIRM", case=False, na=False
        )
    else:
        enriched["is_confirmatory"] = False
    return enriched


def finalize_assay_records(df: pd.DataFrame) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    for column in ["assay_chembl_id", "assay_type", "assay_test_type", "description"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    validated = validate_assays(prepared, context="assay_finalization")
    return validated


PIPELINE_STEPS = [
    StepDefinition("normalize_assay_metadata", normalize_assay_metadata),
    StepDefinition("enrich_assay_flags", enrich_assay_flags),
    StepDefinition("finalize_assay_records", finalize_assay_records),
]


def run_assay_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> RunnerResult:
    """Run the assay postprocessing pipeline."""

    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=ASSAY_SCHEMA,
        pipeline_version=pipeline_version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_STEPS",
    "finalize_assay_records",
    "normalize_assay_metadata",
    "run_assay_pipeline",
    "enrich_assay_flags",
]
