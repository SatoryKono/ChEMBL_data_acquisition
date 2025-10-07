"""Transformation steps for target postprocessing."""
from __future__ import annotations

from collections.abc import Sequence
import re

import pandas as pd
 
from library.postprocess.common import StepDefinition, run_steps
from library.postprocess.common.logging import PipelineRunMetrics

 
from library.pipelines.common.metadata import get_pipeline_version
 
from library.postprocess.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
 
from .schema import TARGET_SCHEMA, validate_targets


def normalize_target_fields(
    df: pd.DataFrame,
    *,
    normalize_taxonomy: bool = False,
    fill_missing_identifiers: bool = False,
) -> pd.DataFrame:
    """Normalize textual fields and optionally harmonise taxonomy metadata."""

    normalized = df.copy(deep=True)

    if "target_chembl_id" in normalized.columns:
        chembl_ids = normalized["target_chembl_id"].astype("string").str.strip()
        normalized["target_chembl_id"] = chembl_ids.str.upper()

        if fill_missing_identifiers:
            missing_mask = chembl_ids.replace({"": pd.NA}).isna()
            fallback_columns = ["chembl_id", "target_id"]
            for column in fallback_columns:
                if not missing_mask.any():
                    break
                if column not in normalized.columns:
                    continue
                fallback = (
                    normalized.loc[missing_mask, column]
                    .astype("string")
                    .str.strip()
                    .replace({"": pd.NA})
                )
                normalized.loc[missing_mask, "target_chembl_id"] = (
                    fallback.str.upper().fillna(normalized.loc[missing_mask, "target_chembl_id"])
                )
                chembl_ids = normalized["target_chembl_id"].astype("string")
                missing_mask = chembl_ids.replace({"": pd.NA}).isna()

    for column in ["pref_name", "organism", "target_type", "target_class", "protein_family", "synonyms"]:
        if column in normalized.columns:
            normalized[column] = (
                normalized[column]
                    .astype("string")
                    .str.strip()
                    .replace({"": pd.NA})
            )

    if normalize_taxonomy:
        taxonomy_columns = [
            "organism",
            "lineage_superkingdom",
            "lineage_phylum",
            "lineage_class",
            "species_group_flag",
        ]
        for column in taxonomy_columns:
            if column in normalized.columns:
                normalized[column] = (
                    normalized[column]
                    .astype("string")
                    .str.strip()
                    .str.replace(r"\s+", " ", regex=True)
                    .replace({"": pd.NA})
                )

        for column in ["taxon_id", "tax_id"]:
            if column in normalized.columns:
                normalized[column] = pd.to_numeric(
                    normalized[column].astype("string").str.strip(), errors="coerce"
                ).astype("Int64")

    return normalized


_SYNONYM_SPLIT_PATTERN = re.compile(r"[|;,]")


def _tokenise_synonyms(value: object) -> list[str]:
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    parts = _SYNONYM_SPLIT_PATTERN.split(text)
    return [part.strip() for part in parts if part and part.strip()]


def enrich_target_synonyms(
    df: pd.DataFrame,
    *,
    synonym_sources: Sequence[str] | None = None,
    preferred_separator: str = ", ",
) -> pd.DataFrame:
    """Aggregate synonyms from configured sources into a deterministic string column."""

    enriched = df.copy(deep=True)
    synonym_columns: list[str] = []

    if "synonyms" in enriched.columns:
        synonym_columns.append("synonyms")

    if synonym_sources:
        for source in synonym_sources:
            column_name = f"{source}_synonyms"
            if column_name in enriched.columns and column_name not in synonym_columns:
                synonym_columns.append(column_name)

    if not synonym_columns:
        return enriched

    def _merge_row(row: pd.Series) -> object:
        seen: set[str] = set()
        ordered_tokens: list[str] = []

        for column in synonym_columns:
            tokens = _tokenise_synonyms(row.get(column))
            for token in tokens:
                marker = token.casefold()
                if marker in seen:
                    continue
                seen.add(marker)
                ordered_tokens.append(token)

        if not ordered_tokens:
            return pd.NA

        return preferred_separator.join(ordered_tokens)

    enriched["synonyms"] = (
        enriched.apply(_merge_row, axis=1)
        .astype("string")
        .replace({"": pd.NA})
    )

    return enriched


def finalize_target_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    sort_by: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`TARGET_SCHEMA`."""

    prepared = df.copy(deep=True)
    for column in TARGET_SCHEMA.required_columns:
        if column not in prepared.columns:
            prepared[column] = pd.Series(pd.NA, index=prepared.index, dtype="string")
    for column in ["target_chembl_id", "pref_name", "target_type"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if enforce_schema:
        validated = validate_targets(prepared, context="target_finalization")
    else:
        ordered_columns: list[str] = []
        if TARGET_SCHEMA.column_order:
            ordered_columns.extend(
                column for column in TARGET_SCHEMA.column_order if column in prepared.columns
            )
        remaining = [
            column for column in prepared.columns if column not in ordered_columns
        ]
        validated = prepared.loc[:, ordered_columns + remaining] if ordered_columns else prepared

    if sort_by:
        sort_columns = [column for column in sort_by if column in validated.columns]
        if sort_columns:
            validated = validated.sort_values(sort_columns, kind="mergesort").reset_index(drop=True)

    return validated


PIPELINE_CONFIG = load_pipeline_config("targets")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def run_target_pipeline(
    df: pd.DataFrame, *, pipeline_version: str | None = None, logger=None
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Run the target postprocessing pipeline and return metrics."""

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
