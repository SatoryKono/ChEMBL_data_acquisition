"""Transformation steps for the test item postprocessing pipeline."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocessing.pipeline.common import run_steps
from library.postprocessing.pipeline.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
from library.postprocessing.pipeline.common.logging import PipelineRunMetrics

from .schema import BOOLEAN_COLUMNS, TESTITEM_SCHEMA, validate_testitems

TRUE_MARKERS = {"true", "1", "yes", "y", "t"}
FALSE_MARKERS = {"false", "0", "no", "n", "f"}


def _normalise_identifier(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.upper()


def _normalise_string(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def _coerce_boolean(series: pd.Series) -> pd.Series:
    def _convert(value: object) -> object:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return pd.NA
        if isinstance(value, bool):
            return value
        text = str(value).strip().lower()
        if not text:
            return pd.NA
        if text in TRUE_MARKERS:
            return True
        if text in FALSE_MARKERS:
            return False
        return pd.NA

    coerced = series.map(_convert)
    return pd.Series(pd.array(coerced, dtype=pd.BooleanDtype()), index=series.index)


def normalize_testitem_records(
    df: pd.DataFrame,
    *,
    uppercase_identifiers: bool = True,
    coerce_booleans: bool = True,
    string_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Normalise identifier casing and canonical string columns."""

    normalized = df.copy(deep=True)

    identifier_columns = (
        "molecule_chembl_id",
        "parent_molecule_chembl_id",
        "salt_chembl_id",
    )
    if uppercase_identifiers:
        for column in identifier_columns:
            if column in normalized.columns:
                normalized[column] = _normalise_identifier(normalized[column])
    else:
        for column in identifier_columns:
            if column in normalized.columns:
                normalized[column] = normalized[column].astype("string").str.strip()

    default_string_columns = (
        "pref_name",
        "molecule_type",
        "structure_type",
        "max_phase",
        "canonical_smiles",
        "standard_inchi",
        "standard_inchi_key",
        "pubchem_canonical_smiles",
        "pubchem_isomeric_smiles",
        "pubchem_inchi",
        "pubchem_inchikey",
        "pubchem_iupac_name",
        "pubchem_molecular_formula",
        "black_box_warning",
        "first_approval",
        "timestamp_utc",
    )

    target_columns = string_columns or default_string_columns
    for column in target_columns:
        if column in normalized.columns:
            normalized[column] = _normalise_string(normalized[column])

    if "pubchem_cid" in normalized.columns:
        normalized["pubchem_cid"] = normalized["pubchem_cid"].astype(object)

    if coerce_booleans:
        for column in BOOLEAN_COLUMNS:
            if column in normalized.columns:
                normalized[column] = _coerce_boolean(normalized[column])

    return normalized


def enrich_testitem_annotations(
    df: pd.DataFrame,
    *,
    fill_missing_columns: bool = True,
) -> pd.DataFrame:
    """Ensure optional annotation columns exist with deterministic defaults."""

    enriched = df.copy(deep=True)
    if not fill_missing_columns:
        return enriched

    row_count = len(enriched.index)
    for column in TESTITEM_SCHEMA.optional_columns:
        if column in enriched.columns:
            continue
        dtype = TESTITEM_SCHEMA.dtypes.get(column)
        if isinstance(dtype, pd.BooleanDtype):
            enriched[column] = pd.Series(
                pd.array([pd.NA] * row_count, dtype=dtype), index=enriched.index
            )
        elif dtype is object:
            enriched[column] = pd.Series(
                [pd.NA] * row_count, index=enriched.index, dtype="object"
            )
        else:
            target_dtype = "string" if dtype is None else dtype
            enriched[column] = pd.Series(
                pd.array([pd.NA] * row_count, dtype=target_dtype), index=enriched.index
            )
    return enriched


def _ensure_column_types(df: pd.DataFrame) -> pd.DataFrame:
    coerced = df.copy(deep=True)
    for column, dtype in TESTITEM_SCHEMA.dtypes.items():
        if column not in coerced.columns:
            continue
        try:
            coerced[column] = coerced[column].astype(dtype)
        except TypeError:
            if dtype is object:
                coerced[column] = coerced[column].astype("object")
            else:
                coerced[column] = coerced[column].astype("string")
    return coerced


def finalize_testitem_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    pipeline_version: str | None = None,
    sort_by: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Validate and reorder the DataFrame according to :data:`TESTITEM_SCHEMA`."""

    prepared = df.copy(deep=True)

    if pipeline_version:
        prepared["pipeline_version"] = prepared.get(
            "pipeline_version",
            pd.Series(
                pd.array([pipeline_version] * len(prepared.index), dtype="string"),
                index=prepared.index,
            ),
        )
        prepared["pipeline_version"] = (
            prepared["pipeline_version"].astype("string").fillna(str(pipeline_version))
        )

    prepared = _ensure_column_types(prepared)

    if enforce_schema:
        validated = validate_testitems(prepared, context="testitems_finalization")
    else:
        validated = prepared

    ordering = sort_by or TESTITEM_SCHEMA.sort_by or ()
    if ordering:
        sort_columns = [column for column in ordering if column in validated.columns]
        if sort_columns:
            validated = validated.sort_values(
                sort_columns, kind="mergesort"
            ).reset_index(drop=True)

    return validated


PIPELINE_CONFIG = load_pipeline_config("testitems")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


def _resolve_pipeline_version(override: str | None) -> str:
    candidate = normalize_pipeline_version(override)
    if candidate is not None:
        return candidate

    config_candidate = normalize_pipeline_version(PIPELINE_CONFIG.pipeline_version)
    if config_candidate is not None:
        return config_candidate

    return get_pipeline_version()


def run_testitem_pipeline(
    df: pd.DataFrame,
    *,
    pipeline_version: str | None = None,
    logger=None,
) -> tuple[pd.DataFrame, PipelineRunMetrics]:
    """Execute the test item postprocessing pipeline and return metrics."""

    resolved_version = _resolve_pipeline_version(pipeline_version)
    return run_steps(
        df,
        PIPELINE_STEPS,
        post_schema=TESTITEM_SCHEMA,
        pipeline_version=resolved_version,
        logger=logger,
    )


__all__ = [
    "PIPELINE_CONFIG",
    "PIPELINE_STEPS",
    "enrich_testitem_annotations",
    "finalize_testitem_records",
    "normalize_testitem_records",
    "run_testitem_pipeline",
]
