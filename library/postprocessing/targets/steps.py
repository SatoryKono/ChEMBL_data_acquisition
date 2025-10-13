"""Transformation steps for target postprocessing."""

from __future__ import annotations

import logging
import re
from collections.abc import Sequence
from functools import lru_cache

import pandas as pd

from library.pipelines.common.metadata import get_pipeline_version
from library.postprocessing.common import run_steps
from library.postprocessing.common.config import (
    load_pipeline_config,
    normalize_pipeline_version,
)
from library.postprocessing.common.logging import PipelineRunMetrics
from library.resources.dictionaries import DictionaryManifestError, get_resource

from .schema import TARGET_SCHEMA, validate_targets

LOGGER = logging.getLogger(__name__)


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
                normalized.loc[
                    missing_mask, "target_chembl_id"
                ] = fallback.str.upper().fillna(
                    normalized.loc[missing_mask, "target_chembl_id"]
                )
                chembl_ids = normalized["target_chembl_id"].astype("string")
                missing_mask = chembl_ids.replace({"": pd.NA}).isna()

    for column in [
        "pref_name",
        "organism",
        "target_type",
        "target_class",
        "protein_family",
        "synonyms",
    ]:
        if column in normalized.columns:
            normalized[column] = (
                normalized[column].astype("string").str.strip().replace({"": pd.NA})
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
    **_: object,
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
        enriched.apply(_merge_row, axis=1).astype("string").replace({"": pd.NA})
    )

    return enriched


@lru_cache(maxsize=1)
def _load_target_types_frame() -> pd.DataFrame:
    """Return cached target classification dictionary."""

    try:
        resource = get_resource("target_types")
    except DictionaryManifestError as exc:
        LOGGER.warning(
            "Target metadata enrichment skipped; resource 'target_types' missing: %s",
            exc,
        )
        return pd.DataFrame()

    try:
        frame = pd.read_csv(resource.path, dtype="string")
    except FileNotFoundError as exc:
        LOGGER.warning(
            "Target metadata enrichment skipped; dictionary file not found: %s",
            exc,
        )
        return pd.DataFrame()
    except ValueError as exc:  # pragma: no cover - defensive guard
        LOGGER.warning(
            "Target metadata enrichment skipped; failed to parse %s: %s",
            resource.path,
            exc,
        )
        return pd.DataFrame()

    required = {"target_chembl_id", "synonyms", "IUPHAR_type", "target_id"}
    missing = required - set(frame.columns)
    if missing:
        LOGGER.warning(
            "Target metadata enrichment skipped; dictionary %s missing columns %s",
            resource.path,
            ", ".join(sorted(missing)),
        )
        return pd.DataFrame()

    normalized = frame.loc[:, ["target_chembl_id", "synonyms", "IUPHAR_type", "target_id"]].copy()
    normalized = normalized.fillna("")
    normalized["target_chembl_id"] = (
        normalized["target_chembl_id"].astype("string").str.strip().str.upper()
    )
    normalized["synonyms"] = normalized["synonyms"].astype("string").str.strip()
    normalized["IUPHAR_type"] = normalized["IUPHAR_type"].astype("string").str.strip()
    normalized["target_id"] = normalized["target_id"].astype("string").str.strip()

    normalized = normalized.rename(
        columns={
            "synonyms": "chembl_synonyms_lookup",
            "IUPHAR_type": "target_class_lookup",
            "target_id": "iuphar_target_id_lookup",
        }
    )
    normalized = normalized.drop_duplicates("target_chembl_id", keep="first")
    return normalized


@lru_cache(maxsize=1)
def _load_iuphar_targets_frame() -> pd.DataFrame:
    """Return cached IUPHAR target dictionary."""

    try:
        resource = get_resource("target_iuphar_target")
    except DictionaryManifestError as exc:
        LOGGER.warning(
            "Target metadata enrichment skipped; resource 'target_iuphar_target' missing: %s",
            exc,
        )
        return pd.DataFrame()

    try:
        frame = pd.read_csv(resource.path, dtype="string")
    except FileNotFoundError as exc:
        LOGGER.warning(
            "Target metadata enrichment skipped; IUPHAR file not found: %s",
            exc,
        )
        return pd.DataFrame()
    except ValueError as exc:  # pragma: no cover - defensive guard
        LOGGER.warning(
            "Target metadata enrichment skipped; failed to parse %s: %s",
            resource.path,
            exc,
        )
        return pd.DataFrame()

    required = {"target_id", "synonyms", "family_name"}
    missing = required - set(frame.columns)
    if missing:
        LOGGER.warning(
            "Target metadata enrichment skipped; IUPHAR dictionary %s missing columns %s",
            resource.path,
            ", ".join(sorted(missing)),
        )
        return pd.DataFrame()

    normalized = frame.loc[:, ["target_id", "synonyms", "family_name"]].copy()
    normalized = normalized.fillna("")
    normalized["target_id"] = normalized["target_id"].astype("string").str.strip()
    normalized["synonyms"] = normalized["synonyms"].astype("string").str.strip()
    normalized["family_name"] = normalized["family_name"].astype("string").str.strip()

    normalized = normalized.rename(
        columns={
            "target_id": "iuphar_target_id_lookup",
            "synonyms": "gtopdb_synonyms_lookup",
            "family_name": "protein_family_lookup",
        }
    )
    normalized = normalized.drop_duplicates("iuphar_target_id_lookup", keep="first")
    return normalized


def _get_target_types_frame() -> pd.DataFrame:
    try:
        return _load_target_types_frame()
    except DictionaryManifestError as exc:  # pragma: no cover - safeguard
        LOGGER.warning("Failed to load target type dictionary: %s", exc)
        return pd.DataFrame()


def _get_iuphar_targets_frame() -> pd.DataFrame:
    try:
        return _load_iuphar_targets_frame()
    except DictionaryManifestError as exc:  # pragma: no cover - safeguard
        LOGGER.warning("Failed to load IUPHAR target dictionary: %s", exc)
        return pd.DataFrame()


def _coalesce_column(frame: pd.DataFrame, target: str, source: str) -> None:
    """Merge ``source`` column into ``target`` preferring existing values."""

    if source not in frame.columns:
        return

    source_series = (
        frame[source]
        .astype("string")
        .str.strip()
        .replace({"": pd.NA})
    )
    if target not in frame.columns:
        frame[target] = source_series
    else:
        frame[target] = (
            frame[target]
            .astype("string")
            .replace({"": pd.NA})
            .combine_first(source_series)
        )
    frame.drop(columns=[source], inplace=True)


def enrich_target_metadata(df: pd.DataFrame, **_: object) -> pd.DataFrame:
    """Enrich targets with classification and synonym metadata."""

    if df.empty:
        return df.copy(deep=True)

    if "target_chembl_id" not in df.columns:
        return df.copy(deep=True)

    types_lookup = _get_target_types_frame()
    if types_lookup.empty:
        return df.copy(deep=True)

    enriched = df.copy(deep=True)
    normalized_ids = (
        enriched["target_chembl_id"].astype("string").str.strip().str.upper()
    )
    enriched["_normalized_target_chembl_id"] = normalized_ids

    joined = enriched.merge(
        types_lookup.rename(columns={"target_chembl_id": "_lookup_target_chembl_id"}),
        left_on="_normalized_target_chembl_id",
        right_on="_lookup_target_chembl_id",
        how="left",
        sort=False,
    )

    iuphar_lookup = _get_iuphar_targets_frame()
    if not iuphar_lookup.empty:
        joined = joined.merge(
            iuphar_lookup,
            on="iuphar_target_id_lookup",
            how="left",
            sort=False,
        )

    _coalesce_column(joined, "chembl_synonyms", "chembl_synonyms_lookup")
    _coalesce_column(joined, "target_class", "target_class_lookup")
    _coalesce_column(joined, "protein_family", "protein_family_lookup")
    _coalesce_column(joined, "gtopdb_synonyms", "gtopdb_synonyms_lookup")

    helper_columns = [
        column
        for column in (
            "_normalized_target_chembl_id",
            "_lookup_target_chembl_id",
            "iuphar_target_id_lookup",
        )
        if column in joined.columns
    ]
    if helper_columns:
        joined = joined.drop(columns=helper_columns)

    return joined


def finalize_target_records(
    df: pd.DataFrame,
    *,
    enforce_schema: bool = True,
    sort_by: Sequence[str] | None = None,
    **_: object,
) -> pd.DataFrame:
    """Validate and order the DataFrame according to :data:`TARGET_SCHEMA`."""

    prepared = df.copy(deep=True)

    # Normalise a few legacy aliases that may appear in historical exports.
    if "target_type" not in prepared.columns:
        fallback_columns = (
            "relationship",
            "target_type_description",
            "target_type_desc",
            "target_type_name",
            "targettype",
            "type",
        )
        for alias in fallback_columns:
            if alias in prepared.columns:
                prepared["target_type"] = prepared[alias]
                break

    required_columns = set(TARGET_SCHEMA.required_columns)
    optional_columns = set(TARGET_SCHEMA.optional_columns)
    expected_string_columns = required_columns.union(
        {column for column, dtype in TARGET_SCHEMA.dtypes.items() if dtype == "string"}
    )

    for column in sorted(required_columns | optional_columns):
        if column not in prepared.columns:
            prepared[column] = pd.Series(pd.NA, index=prepared.index, dtype="string")

    for column in expected_string_columns:
        if column in prepared.columns:
            prepared[column] = prepared[column].astype("string").replace({"": pd.NA})

    prepared = _populate_target_type(prepared)

    if enforce_schema:
        validated = validate_targets(prepared, context="target_finalization")
    else:
        ordered_columns: list[str] = []
        if TARGET_SCHEMA.column_order:
            ordered_columns.extend(
                column
                for column in TARGET_SCHEMA.column_order
                if column in prepared.columns
            )
        remaining = [
            column for column in prepared.columns if column not in ordered_columns
        ]
        validated = (
            prepared.loc[:, ordered_columns + remaining]
            if ordered_columns
            else prepared
        )

    if sort_by:
        sort_columns = [column for column in sort_by if column in validated.columns]
        if sort_columns:
            validated = validated.sort_values(
                sort_columns, kind="mergesort"
            ).reset_index(drop=True)

    return validated


PIPELINE_CONFIG = load_pipeline_config("targets")
PIPELINE_STEPS = PIPELINE_CONFIG.step_definitions()


@lru_cache(maxsize=1)
def _load_target_type_lookup() -> pd.Series:
    """Return cached mapping of target identifiers to target type labels."""

    resource = get_resource("target_types")
    try:
        frame = pd.read_csv(
            resource.path,
            usecols=["target_chembl_id", "type"],
            dtype="string",
        )
    except ValueError as exc:  # pragma: no cover - defensive guard
        LOGGER.warning(
            "Failed to load target type dictionary (%s): %s",
            resource.path,
            exc,
        )
        return pd.Series(dtype="string")

    if frame.empty or "target_chembl_id" not in frame or "type" not in frame:
        LOGGER.warning(
            "Target type dictionary %s is missing required columns",
            resource.path,
        )
        return pd.Series(dtype="string")

    normalised = frame.copy()
    normalised["target_chembl_id"] = (
        normalised["target_chembl_id"]
        .astype("string")
        .str.strip()
        .str.upper()
        .replace({"": pd.NA})
    )
    normalised["type"] = (
        normalised["type"].astype("string").str.strip().replace({"": pd.NA})
    )

    normalised = normalised.dropna(subset=["target_chembl_id"])
    if normalised.empty:
        return pd.Series(dtype="string")

    normalised = normalised.drop_duplicates("target_chembl_id", keep="first")
    lookup = normalised.set_index("target_chembl_id")["type"].astype("string")
    return lookup


def _get_target_type_lookup() -> pd.Series:
    try:
        return _load_target_type_lookup()
    except (FileNotFoundError, DictionaryManifestError) as exc:
        LOGGER.warning(
            "Target type dictionary resource 'target_types' could not be loaded: %s",
            exc,
        )
        return pd.Series(dtype="string")


def _populate_target_type(df: pd.DataFrame) -> pd.DataFrame:
    """Fill missing ``target_type`` values using the bundled dictionary."""

    if "target_type" not in df.columns or "target_chembl_id" not in df.columns:
        return df

    missing_mask = df["target_type"].isna()
    if not missing_mask.any():
        return df

    lookup = _get_target_type_lookup()
    if lookup.empty:
        return df

    identifiers = (
        df["target_chembl_id"]
        .astype("string")
        .str.strip()
        .str.upper()
        .replace({"": pd.NA})
    )
    mapped = identifiers[missing_mask].map(lookup)
    if mapped.isna().all():
        return df

    df.loc[missing_mask, "target_type"] = mapped.astype("string").replace({"": pd.NA})
    return df


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
    "enrich_target_metadata",
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "enrich_target_synonyms",
]
