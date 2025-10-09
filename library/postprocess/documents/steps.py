"""Transformation steps for document postprocessing."""
from __future__ import annotations

import unicodedata

import pandas as pd
 
from library.postprocess.common import run_steps
from library.postprocess.common.logging import PipelineRunMetrics
 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 

from .schema import DOCUMENT_SCHEMA, validate_documents


def normalize_document_fields(
    df: pd.DataFrame,
    *,
    trim_whitespace: bool = True,
    normalise_unicode: bool = False,
) -> pd.DataFrame:
    """Normalize whitespace and casing for textual document fields."""

    normalized = df.copy(deep=True)
    normalized.columns = [col.strip().lower() for col in normalized.columns]

    textual_columns = ["title", "journal", "doc_type"]
    for column in textual_columns:
        if column not in normalized.columns:
            continue
        series = normalized[column].astype("string")
        if trim_whitespace:
            series = series.str.strip().str.replace(r"\s+", " ", regex=True)
        if normalise_unicode:
            series = series.map(
                lambda value: unicodedata.normalize("NFKC", value)
                if not pd.isna(value)
                else value
            )
        normalized[column] = series.replace({"": pd.NA})

    return normalized


_YEAR_VALID_RANGE = (1600, 2100)


def _coerce_year(series: pd.Series) -> pd.Series:
    """Return a nullable integer series containing plausible year values."""

    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.empty:
        return pd.Series(pd.NA, index=series.index, dtype="Int64")

    lower, upper = _YEAR_VALID_RANGE
    mask = numeric.between(lower, upper)
    cleaned = numeric.where(mask)
    return cleaned.round().astype("Int64")


def enrich_document_publication_year(
    df: pd.DataFrame,
    *,
    fallback_year: int | None = None,
    prefer_doi_year: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Create a deterministic ``publication_year`` column.

    Parameters
    ----------
    df:
        Input DataFrame produced by earlier pipeline steps.
    fallback_year:
        Optional numeric year to use when no source provides a valid value.
    prefer_doi_year:
        When ``True``, favour years originating from external DOI metadata over the
        ChEMBL export column.
    """

    if fallback_year is not None:
        try:
            fallback_value = int(fallback_year)
        except (TypeError, ValueError) as exc:
            raise ValueError("fallback_year must be an integer") from exc
        if not (_YEAR_VALID_RANGE[0] <= fallback_value <= _YEAR_VALID_RANGE[1]):
            raise ValueError(
                "fallback_year must be within the supported range "
                f"{_YEAR_VALID_RANGE[0]}-{_YEAR_VALID_RANGE[1]}"
            )
    else:
        fallback_value = None

    enriched = df.copy(deep=True)
    # Determine priority order for candidate columns.
    supplemental_sources = [
        "crossref.year",
        "crossref.published",
        "openalex.publication_year",
        "openalex.year",
        "scholar.year",
        "pubmed.yearcompleted",
        "pubmed.yearrevised",
        "pubmed.year",
    ]
    primary_sources = ["year"]

    if prefer_doi_year:
        ordered_candidates = supplemental_sources + primary_sources
    else:
        ordered_candidates = primary_sources + supplemental_sources

    publication_year = pd.Series(pd.NA, index=enriched.index, dtype="Int64")

    for column in ordered_candidates:
        if column not in enriched.columns:
            continue
        candidate = _coerce_year(enriched[column])
        if candidate.isna().all():
            continue
        publication_year = publication_year.fillna(candidate)
        if not publication_year.isna().any():
            break

    if fallback_value is not None:
        fallback_series = pd.Series(fallback_value, index=enriched.index, dtype="Int64")
        publication_year = publication_year.fillna(fallback_series)

    enriched["publication_year"] = publication_year
    return enriched


def finalize_document_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    ensure_unique_ids: bool = False,
    **_: object,
) -> pd.DataFrame:
    """Normalise terminal document rows and optionally validate the schema.

    Parameters
    ----------
    df:
        Input frame containing the document rows produced by the preceding
        enrichment steps.
    enforce_schema:
        When ``True`` (the default) the output is validated against
        :data:`DOCUMENT_SCHEMA`. Setting this flag to ``False`` can be useful
        in exploratory runs where partial data is expected.
    ensure_unique_ids:
        Drop duplicated ``document_chembl_id`` values while keeping the first
        occurrence. The operation is deterministic to guarantee idempotent
        pipeline re-runs.
    """

    prepared = df.copy(deep=True)

    identifier_column = "document_chembl_id"
    id_fallback_columns = (
        identifier_column,
        "chembl.document_chembl_id",
    )

    if identifier_column in prepared.columns:
        document_ids = prepared[identifier_column].astype("string")
    else:
        document_ids = pd.Series(pd.NA, index=prepared.index, dtype="string")

    document_ids = document_ids.replace({"": pd.NA})

    for fallback in id_fallback_columns:
        if fallback == identifier_column or fallback not in prepared.columns:
            continue
        candidate = prepared[fallback].astype("string").replace({"": pd.NA})
        document_ids = document_ids.fillna(candidate)

    prepared[identifier_column] = document_ids

    required_string_columns = [identifier_column, "title", "doc_type"]
    for column in required_string_columns:
        if column not in prepared.columns:
            prepared[column] = pd.Series(pd.NA, index=prepared.index, dtype="string")
        else:
            prepared[column] = prepared[column].astype("string").replace({"": pd.NA})

    if "publication_year" in prepared.columns:
        prepared["publication_year"] = prepared["publication_year"].astype("Int64")

    if ensure_unique_ids and "document_chembl_id" in prepared.columns:
        prepared = (
            prepared.drop_duplicates(subset=["document_chembl_id"], keep="first")
            .reset_index(drop=True)
        )
    else:
        prepared = prepared.reset_index(drop=True)

    if not enforce_schema:
        ordered = prepared
        column_order = DOCUMENT_SCHEMA.column_order or ()
        if column_order:
            ordered_columns = [col for col in column_order if col in ordered.columns]
            remaining = [col for col in ordered.columns if col not in ordered_columns]
            if ordered_columns:
                ordered = ordered.loc[:, ordered_columns + remaining]

        sort_by = DOCUMENT_SCHEMA.sort_by or ()
        if sort_by:
            sort_columns = [col for col in sort_by if col in ordered.columns]
            if sort_columns:
                ordered = ordered.sort_values(sort_columns, kind="mergesort").reset_index(
                    drop=True
                )
        return ordered

    validated = validate_documents(prepared, context="document_finalization")
    return validated


PIPELINE_CONFIG = load_pipeline_config("documents")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_document_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the document postprocessing pipeline and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=DOCUMENT_SCHEMA,
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
    "finalize_document_records",
    "normalize_document_fields",
    "run_document_pipeline",
    "enrich_document_publication_year",
]
