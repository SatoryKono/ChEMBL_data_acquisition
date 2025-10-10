"""Utilities for preparing target exports for schema validation."""

from __future__ import annotations

from typing import Iterable

import pandas as pd

from library.schemas import TargetsSchema
from library.schemas.targets import TARGETS_COLUMN_ORDER

TARGETS_REQUIRED_COLUMNS: set[str] = {
    name for name, column in TargetsSchema.columns.items() if column.required
}

TARGETS_OPTIONAL_COLUMNS: list[str] = [
    column for column in TARGETS_COLUMN_ORDER if column not in TARGETS_REQUIRED_COLUMNS
]

TARGETS_OBJECT_COLUMNS: set[str] = {
    name
    for name, column in TargetsSchema.columns.items()
    if str(column.dtype) == "object"
}


def _placeholder_series(index: pd.Index, *, value: str = "-") -> pd.Series:
    if not len(index):
        return pd.Series(dtype=object)
    return pd.Series([value] * len(index), index=index, dtype=object)


def _coerce_optional_columns(
    df: pd.DataFrame, missing_optional: Iterable[str]
) -> dict[str, pd.Series]:
    placeholders: dict[str, pd.Series] = {}
    for column in missing_optional:
        placeholders[column] = _placeholder_series(df.index)
    return placeholders


def prepare_targets_for_schema(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Align *df* to :data:`TargetsSchema` and report missing columns."""

    present_columns = set(df.columns)
    missing_required = TARGETS_REQUIRED_COLUMNS - present_columns
    missing_optional = {
        column for column in TARGETS_OPTIONAL_COLUMNS if column not in present_columns
    }

    prepared = df.copy()

    if missing_optional:
        prepared = prepared.assign(
            **_coerce_optional_columns(prepared, sorted(missing_optional))
        )

    extras: dict[str, pd.Series] = {}

    if "uniprot_id_primary" not in prepared.columns and "uniprot_id" in prepared.columns:
        extras["uniprot_id_primary"] = prepared["uniprot_id"].astype(object)

    if "mapping_uniprot_id" in prepared.columns:
        extras["mapping_uniprot_id"] = prepared["mapping_uniprot_id"].astype(object)

    ordered_columns = list(TARGETS_COLUMN_ORDER)
    for column in extras:
        if column not in ordered_columns:
            ordered_columns.append(column)

    prepared = prepared.reindex(columns=ordered_columns)

    if extras:
        for column, series in extras.items():
            prepared[column] = series.reindex(prepared.index)

    for column in TARGETS_OBJECT_COLUMNS & set(prepared.columns):
        prepared[column] = prepared[column].astype(object)

    return prepared, missing_required, set(missing_optional)


__all__ = [
    "prepare_targets_for_schema",
    "TARGETS_OPTIONAL_COLUMNS",
    "TARGETS_REQUIRED_COLUMNS",
    "TARGETS_OBJECT_COLUMNS",
]
