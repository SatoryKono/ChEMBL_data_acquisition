"""Transformation steps for target postprocessing."""
from __future__ import annotations

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


def enrich_target_synonyms(
    df: pd.DataFrame,
    *,
    preferred_separator: str = ", ",
    synonym_sources: tuple[str, ...] | list[str] | None = None,
) -> pd.DataFrame:
    """Ensure the ``synonyms`` column is deterministically ordered.

    Parameters
    ----------
    df:
        Input frame.
    preferred_separator:
        Separator used to join the normalised synonym tokens.  Defaults to
        ``", "`` for backwards compatibility with the legacy behaviour.
    synonym_sources:
        A declarative list of sources in priority order.  At present the
        pipeline delivers a single denormalised ``synonyms`` column.  The
        parameter is therefore accepted for compatibility with the pipeline
        configuration, but no source specific handling is required yet.  The
        value is ignored other than for truthiness checks to keep the
        signature aligned with the YAML configuration.
    """

    enriched = df.copy(deep=True)
    if "synonyms" in enriched.columns:
        tokens = (
            enriched["synonyms"]
            .fillna("")
            .astype("string")
            .str.split(",")
            .apply(
                lambda values: sorted(
                    {part.strip() for part in values if isinstance(part, str) and part.strip()}
                )
            )
        )
        enriched["synonyms"] = tokens.apply(
            lambda values: preferred_separator.join(values)
        )
    # ``synonym_sources`` is currently unused; the argument is accepted so the
    # pipeline configuration does not trigger warnings about unsupported
    # parameters.  A simple truthiness check makes the intent explicit and keeps
    # static analysers quiet without altering behaviour.
    if synonym_sources:
        _ = tuple(synonym_sources)
    return enriched


def finalize_target_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    sort_by: tuple[str, ...] | list[str] | None = None,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`TARGET_SCHEMA`.

    The pipeline configuration historically supplied ``enforce_schema`` and
    ``sort_by`` parameters which were silently ignored.  Accepting the
    arguments here keeps the step behaviour aligned with the declarative YAML
    configuration and prevents noisy warnings during pipeline execution.
    """

    prepared = df.copy(deep=True)

    # Ensure all required columns are present before performing any validation.
    for column in TARGET_SCHEMA.required_columns:
        if column not in prepared.columns:
            dtype = TARGET_SCHEMA.dtypes.get(column, "string")
            prepared[column] = pd.Series(pd.NA, index=prepared.index, dtype=dtype)

    # Coerce the canonical string columns to pandas' ``string`` dtype so that
    # downstream schema validation observes consistent types across pandas
    # releases.
    for column in ["target_chembl_id", "pref_name", "target_type"]:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string")

    if enforce_schema:
        validated = validate_targets(prepared, context="target_finalization")
    else:
        validated = prepared.copy(deep=True)
        column_order = TARGET_SCHEMA.column_order
        if column_order:
            ordered_columns = [
                column for column in column_order if column in validated.columns
            ]
            remaining = [
                column
                for column in validated.columns
                if column not in ordered_columns
            ]
            validated = validated.loc[:, ordered_columns + remaining]

    sort_columns: list[str] = []
    if sort_by:
        sort_columns = [column for column in sort_by if column in validated.columns]
    elif TARGET_SCHEMA.sort_by:
        sort_columns = [
            column for column in TARGET_SCHEMA.sort_by if column in validated.columns
        ]
    if sort_columns:
        validated = validated.sort_values(sort_columns, kind="mergesort").reset_index(
            drop=True
        )

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
