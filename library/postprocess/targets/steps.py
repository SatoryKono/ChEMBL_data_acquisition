"""Transformation steps for target postprocessing."""
from __future__ import annotations

from typing import Sequence

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.common import PipelineRunMetadata, run_steps
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)

from .schema import TARGET_SCHEMA, validate_targets


def normalize_target_fields(
    df: pd.DataFrame,
    *,
    normalize_taxonomy: bool = True,
    fill_missing_identifiers: bool = True,
    **_: object,
) -> pd.DataFrame:
    """Normalize leading/trailing whitespace and ensure upper-case identifiers."""

    normalized = df.copy(deep=True)

    if fill_missing_identifiers and "target_chembl_id" not in normalized.columns:
        normalized["target_chembl_id"] = pd.Series(dtype="string")

    if "target_chembl_id" in normalized.columns:
        normalized["target_chembl_id"] = (
            normalized["target_chembl_id"].astype("string").str.strip().str.upper()
        )

    for column in ["pref_name", "organism"]:
        if column not in normalized.columns:
            continue
        series = normalized[column].astype("string")
        if normalize_taxonomy:
            series = series.str.strip()
        normalized[column] = series

    return normalized


def enrich_target_synonyms(
    df: pd.DataFrame,
    *,
    synonym_sources: Sequence[str] | None = None,
    preferred_separator: str = ", ",
    **_: object,
) -> pd.DataFrame:
    """Ensure synonyms column is deterministically ordered."""

    enriched = df.copy(deep=True)
    if "synonyms" in enriched.columns:
        enriched["synonyms"] = (
            enriched["synonyms"]
            .fillna("")
            .astype("string")
            .apply(
                lambda value: preferred_separator.join(
                    sorted(
                        part.strip()
                        for part in value.split(",")
                        if part.strip()
                    )
                )
            )
        )
    return enriched


def finalize_target_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    sort_by: Sequence[str] | None = None,
    **_: object,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`TARGET_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in ["target_chembl_id", "pref_name", "target_type"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if not enforce_schema:
        result = prepared
    else:
        result = validate_targets(prepared, context="target_finalization")

    if sort_by:
        available = [column for column in sort_by if column in result.columns]
        if available:
            result = result.sort_values(available, kind="mergesort").reset_index(
                drop=True
            )

    return result


PIPELINE_CONFIG = load_pipeline_config("targets")
PIPELINE_STEPS = PIPELINE_CONFIG.steps


def run_target_pipeline(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetadata]:
    """Run the target postprocessing pipeline."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=TARGET_SCHEMA,
        pipeline_version=resolved_version,
        logger=logger,
    )


def _resolve_pipeline_version(override: str | None) -> str:
    candidate = normalize_pipeline_version(override)
    if candidate is not None:
        return candidate

    config_candidate = normalize_pipeline_version(PIPELINE_CONFIG.pipeline_version)
    if config_candidate is not None:
        return config_candidate

    return get_pipeline_version()


__all__ = [
    "PIPELINE_CONFIG",
    "PIPELINE_STEPS",
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "enrich_target_synonyms",
]
