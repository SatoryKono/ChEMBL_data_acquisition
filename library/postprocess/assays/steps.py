"""Transformation steps for assay postprocessing."""
from __future__ import annotations

import re
import logging
from typing import Sequence

import pandas as pd

 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import ASSAY_SCHEMA, validate_assays


_LOGGER = logging.getLogger(__name__)


def normalize_assay_metadata(
    df: pd.DataFrame,
    *,
    uppercase_categories: bool = True,
    strip_whitespace: bool = True,
    collapse_internal_whitespace: bool | None = None,
    **extra_kwargs: object,
) -> pd.DataFrame:
    """Normalize string-based assay descriptors.

    Parameters
    ----------
    df:
        Input frame to normalize.
    uppercase_categories:
        When ``True`` (the default) categorical assay descriptors are converted to
        upper-case.  Set to ``False`` to preserve the original casing while still
        applying other configured clean-up rules.
    strip_whitespace:
        Remove leading and trailing whitespace around categorical values.  Enabled
        by default for backwards compatibility with historical exports.
    collapse_internal_whitespace:
        Normalise consecutive whitespace characters inside categorical values.  If
        ``None`` the value defaults to ``strip_whitespace`` to preserve the legacy
        behaviour where both operations happened together.  Pass ``True`` or
        ``False`` explicitly to override the coupling.
    """

    if extra_kwargs:
        _LOGGER.debug(
            "normalize_assay_metadata ignoring unsupported parameters: %s",
            ", ".join(sorted(extra_kwargs)),
        )

    normalized = df.copy(deep=True)
    if collapse_internal_whitespace is None:
        collapse_internal_whitespace = strip_whitespace

    for column in ["assay_type", "assay_test_type", "assay_format"]:
        if column not in normalized.columns:
            continue

        series = normalized[column].astype("string")
        if strip_whitespace:
            series = series.str.strip()
        if collapse_internal_whitespace:
            series = series.str.replace(r"\s+", " ", regex=True)
        if uppercase_categories:
            series = series.str.lower()
        normalized[column] = series
    return normalized


def _prepare_confirmatory_terms(
    terms: Sequence[str] | None,
) -> tuple[str, ...]:
    """Return cleaned confirmatory term patterns."""

    if not terms:
        return ()

    prepared: list[str] = []
    for term in terms:
        text = str(term).strip()
        if text:
            prepared.append(text)
    return tuple(prepared)


def enrich_assay_flags(
    df: pd.DataFrame,
    *,
    confirmatory_terms: Sequence[str] | None = None,
    default_flag: bool = False,
) -> pd.DataFrame:
    """Introduce confirmatory flag based on assay type information."""

    enriched = df.copy(deep=True)
    prepared_terms = _prepare_confirmatory_terms(confirmatory_terms)
    base_flag = pd.Series(bool(default_flag), index=enriched.index, dtype="bool")

    type_series = enriched.get("assay_type")
    if type_series is None:
        enriched["is_confirmatory"] = base_flag
        return enriched

    category = type_series.fillna("").astype("string")
    if prepared_terms:
        pattern = "|".join(re.escape(term) for term in prepared_terms)
        matches = category.str.contains(pattern, case=False, regex=True)
        enriched["is_confirmatory"] = matches.fillna(bool(default_flag)).astype(bool)
    else:
        enriched["is_confirmatory"] = base_flag
    return enriched


def finalize_assay_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    normalize_identifiers: bool = True,
    identifier_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Apply schema validation and deterministic ordering."""

    prepared = df.copy(deep=True)
    if identifier_columns is None:
        identifier_columns = ("assay_chembl_id", "target_chembl_id")

    if normalize_identifiers:
        for column in identifier_columns:
            if column in prepared.columns:
                prepared[column] = (
                    prepared[column]
                    .astype("string")
                    .str.strip()
                    .str.upper()
                )

    for column in ["assay_chembl_id", "assay_type", "assay_test_type", "description"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if not enforce_schema:
        return prepared

    validated = validate_assays(prepared, context="assay_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("assays")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_assay_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the assay postprocessing pipeline and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=ASSAY_SCHEMA,
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
    "finalize_assay_records",
    "normalize_assay_metadata",
    "run_assay_pipeline",
    "enrich_assay_flags",
]
