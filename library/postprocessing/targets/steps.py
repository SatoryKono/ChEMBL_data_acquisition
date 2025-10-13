"""Transformation steps for target postprocessing."""

from __future__ import annotations

import json
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


def _normalise_scalar(value: object) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass
    text = str(value).strip()
    if not text:
        return ""
    lowered = text.casefold()
    if lowered in {"nan", "none", "null", "-"}:
        return ""
    return text


def _parse_protein_classifications(value: object) -> list[tuple[int | None, str]]:
    if value is None:
        return []
    if isinstance(value, str):
        text = value.strip()
        if not text or text.casefold() in {"nan", "none", "null", "-"}:
            return []
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError:
            tokens = [token.strip() for token in re.split(r"[|>;]", text) if token.strip()]
            return [
                (index + 1, token)
                for index, token in enumerate(tokens)
            ]
    elif isinstance(value, (list, tuple)):
        parsed = list(value)
    elif isinstance(value, dict):
        parsed = [value]
    else:
        return []

    records: list[tuple[int | None, str]] = []

    for item in parsed:
        if not isinstance(item, dict):
            continue

        candidate = item.get("protein_classification")
        level: int | None = None

        if isinstance(candidate, dict):
            name = (
                candidate.get("pref_name")
                or candidate.get("protein_classification")
                or candidate.get("short_name")
                or candidate.get("name")
            )
            level_value = (
                candidate.get("class_level")
                or candidate.get("protein_classification_level")
                or candidate.get("level")
            )
        else:
            name = candidate
            level_value = None

        if not name:
            name = (
                item.get("pref_name")
                or item.get("protein_classification")
                or item.get("classification")
                or item.get("name")
            )
        if level_value is None:
            level_value = (
                item.get("class_level")
                or item.get("protein_classification_level")
                or item.get("level")
            )

        if name is None:
            continue

        try:
            level = int(level_value) if level_value is not None else None
        except (TypeError, ValueError):
            level = None

        normalised_name = _normalise_scalar(name)
        if not normalised_name:
            continue

        records.append((level, normalised_name))

    return records


def _select_classification(
    entries: list[tuple[int | None, str]],
    preferred_levels: Sequence[int],
) -> str:
    by_level: dict[int | None, list[str]] = {}
    by_level_seen: dict[int | None, set[str]] = {}
    for level, name in entries:
        if not name:
            continue
        key = level if level in preferred_levels else level
        values = by_level.setdefault(key, [])
        seen = by_level_seen.setdefault(key, set())
        marker = name.casefold()
        if marker in seen:
            continue
        seen.add(marker)
        values.append(name)

    for level in preferred_levels:
        values = by_level.get(level)
        if values:
            return "; ".join(values)

    # fall back to any available level preserving discovery order
    seen: set[str] = set()
    ordered: list[str] = []
    for _, name in entries:
        marker = name.casefold()
        if marker in seen:
            continue
        seen.add(marker)
        ordered.append(name)
    return "; ".join(ordered)


def _derive_classification_labels(row: pd.Series) -> tuple[str, str]:
    entries = _parse_protein_classifications(row.get("protein_classifications"))

    target_class = _select_classification(entries, preferred_levels=(1,)) if entries else ""
    protein_family = _select_classification(entries, preferred_levels=(2, 3)) if entries else ""

    if not target_class:
        for column in ("protein_class_pred_L1", "IUPHAR_class", "target_type"):
            candidate = _normalise_scalar(row.get(column))
            if candidate:
                target_class = candidate
                break

    if not protein_family:
        for column in ("protein_class_pred_L2", "IUPHAR_subclass", "IUPHAR_type"):
            candidate = _normalise_scalar(row.get(column))
            if candidate:
                protein_family = candidate
                break

    return target_class, protein_family


def _merge_synonym_tokens(row: pd.Series, columns: Sequence[str]) -> str:
    seen: set[str] = set()
    tokens: list[str] = []
    for column in columns:
        value = row.get(column)
        for token in _tokenise_synonyms(value):
            marker = token.casefold()
            if marker in seen:
                continue
            seen.add(marker)
            tokens.append(token)
    return "|".join(tokens)


def _fill_synonym_column(
    df: pd.DataFrame, column: str, sources: Sequence[str]
) -> None:
    available = [source for source in sources if source in df.columns]
    if not available:
        if column not in df.columns:
            df[column] = pd.Series(pd.NA, index=df.index, dtype="string")
        return

    derived = (
        df.apply(lambda row: _merge_synonym_tokens(row, available), axis=1)
        .astype("string")
        .replace({"": pd.NA})
    )

    if column in df.columns:
        current = df[column].astype("string").replace({"": pd.NA})
    else:
        current = pd.Series(pd.NA, index=df.index, dtype="string")

    df[column] = current.fillna(derived)


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

    classification_sources = {
        "protein_classifications",
        "protein_class_pred_L1",
        "protein_class_pred_L2",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_type",
        "target_type",
    }
    if classification_sources.intersection(normalized.columns):
        derived = normalized.apply(
            lambda row: _derive_classification_labels(row), axis=1, result_type="expand"
        )
        derived.columns = ["_derived_target_class", "_derived_protein_family"]

        if "target_class" in normalized.columns:
            current_class = (
                normalized["target_class"].astype("string").replace({"": pd.NA})
            )
        else:
            current_class = pd.Series(pd.NA, index=normalized.index, dtype="string")

        if "protein_family" in normalized.columns:
            current_family = (
                normalized["protein_family"].astype("string").replace({"": pd.NA})
            )
        else:
            current_family = pd.Series(pd.NA, index=normalized.index, dtype="string")

        normalized["target_class"] = current_class.fillna(
            derived["_derived_target_class"].astype("string").replace({"": pd.NA})
        )
        normalized["protein_family"] = current_family.fillna(
            derived["_derived_protein_family"].astype("string").replace({"": pd.NA})
        )

        normalized = normalized.drop(
            columns=["_derived_target_class", "_derived_protein_family"], errors="ignore"
        )

    _fill_synonym_column(
        normalized,
        "chembl_synonyms",
        [
            "chembl_synonyms",
            "protein_synonym_list",
            "protein_name_canonical",
            "protein_name_alt",
            "pref_name",
            "gene_symbol",
            "gene_symbol_list",
        ],
    )
    _fill_synonym_column(
        normalized,
        "gtopdb_synonyms",
        ["gtopdb_synonyms", "gtop_synonyms"],
    )
    _fill_synonym_column(
        normalized,
        "synonyms",
        [
            "synonyms",
            "protein_synonym_list",
            "protein_name_canonical",
            "protein_name_alt",
            "pref_name",
            "gene_symbol",
            "gene_symbol_list",
        ],
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
    "finalize_target_records",
    "normalize_target_fields",
    "run_target_pipeline",
    "enrich_target_synonyms",
]
