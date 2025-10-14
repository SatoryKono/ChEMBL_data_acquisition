"""Fetch, normalise and validate activity data from the ChEMBL API."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final

import pandas as pd

from library.clients.chembl_client import ChemblClient
from library.schemas.activity_schema import validate_activity
from library.utils.logging import get_logger
from library.utils.qc_report import TableQualityProfiler, build_reports_from_profiler

logger = get_logger(__name__)

ACTIVITY_COLUMN_ORDER: Final[tuple[str, ...]] = (
    "activity_id",
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
)

_API_FIELDS: Final[tuple[str, ...]] = (
    "activity_id",
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
)


@dataclass(frozen=True, slots=True)
class ActivityData:
    """Container bundling the dataset and derived QA artefacts."""

    dataset: pd.DataFrame
    quality_report: pd.DataFrame
    correlation_report: pd.DataFrame
    qc_summary: dict[str, Any]


def fetch_normalize_activity(
    limit: int,
    *,
    chembl_client: ChemblClient | None = None,
) -> pd.DataFrame:
    """Return the validated activity dataset fetched from ChEMBL."""

    if limit <= 0:
        raise ValueError("limit must be positive")

    client = chembl_client if chembl_client is not None else ChemblClient()
    created_client = chembl_client is None
    try:
        records = _fetch_activity_payload(client, limit)
    finally:
        if created_client:
            close = getattr(client, "close", None)
            if callable(close):
                close()

    frame = pd.DataFrame(records, columns=ACTIVITY_COLUMN_ORDER)
    normalised = normalize_activity_frame(frame)
    validated = validate_activity(normalised)
    return validated


def normalize_activity_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` coerced to the canonical schema and ordering."""

    working = frame.copy(deep=True)
    for column in ACTIVITY_COLUMN_ORDER:
        if column not in working.columns:
            working[column] = pd.Series([pd.NA] * len(working))

    working = working.loc[:, ACTIVITY_COLUMN_ORDER]

    working["activity_id"] = pd.Series(
        pd.array(pd.to_numeric(working["activity_id"], errors="coerce"), dtype="Int64"),
        index=working.index,
    )

    for name in ("assay_chembl_id", "target_chembl_id", "standard_type"):
        working[name] = _normalise_string_column(working[name])

    working["standard_relation"] = _normalise_string_column(
        working["standard_relation"], uppercase=True
    )

    working["standard_units"] = _normalise_string_column(
        working["standard_units"], uppercase=True
    )

    numeric = pd.to_numeric(working["standard_value"], errors="coerce")
    working["standard_value"] = pd.Series(
        pd.array(numeric, dtype="Float64"), index=working.index
    )

    working = working.dropna(subset=["activity_id"])
    working = working.drop_duplicates(subset=["activity_id"]).sort_values(
        by=["activity_id", "assay_chembl_id", "target_chembl_id"],
        kind="mergesort",
    )
    working = working.reset_index(drop=True)
    return working


def generate_activity_reports(frame: pd.DataFrame) -> ActivityData:
    """Return the dataset together with QC and correlation artefacts."""

    profiler = TableQualityProfiler()
    profiler.consume(frame)
    quality_report, correlation_report = build_reports_from_profiler(profiler)
    qc_summary = {
        "row_count": int(frame.shape[0]),
        "column_count": int(frame.shape[1]),
        "non_null_ratio": (
            float(frame.notna().mean().mean()) if not frame.empty else float("nan")
        ),
    }
    if math.isnan(qc_summary["non_null_ratio"]):
        qc_summary["non_null_ratio"] = 0.0
    return ActivityData(
        dataset=frame,
        quality_report=quality_report,
        correlation_report=correlation_report,
        qc_summary=qc_summary,
    )


def _fetch_activity_payload(client: ChemblClient, limit: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    offset = 0
    while len(records) < limit:
        page_limit = min(200, limit - len(records))
        payload = client.list_activities(
            limit=page_limit, offset=offset, fields=_API_FIELDS
        )
        items = _extract_activity_items(payload)
        if not items:
            break
        records.extend(_prepare_activity_rows(items))
        offset += len(items)
        if len(records) >= limit:
            break
        page_meta = payload.get("page_meta") if isinstance(payload, Mapping) else None
        next_url = page_meta.get("next") if isinstance(page_meta, Mapping) else None
        if not next_url:
            break
    return records[:limit]


def _extract_activity_items(payload: Mapping[str, Any]) -> Sequence[Mapping[str, Any]]:
    for key in ("activities", "items", "data", "results"):
        value = payload.get(key)
        if isinstance(value, Sequence):
            return [item for item in value if isinstance(item, Mapping)]
    return []


def _prepare_activity_rows(items: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for entry in items:
        activity_id = entry.get("activity_id") or entry.get("id")
        if activity_id is None:
            continue
        activity_token = str(activity_id).strip()
        if not activity_token:
            continue
        rows.append(
            {
                "activity_id": activity_token,
                "assay_chembl_id": entry.get("assay_chembl_id"),
                "target_chembl_id": entry.get("target_chembl_id"),
                "standard_type": entry.get("standard_type"),
                "standard_relation": entry.get("standard_relation"),
                "standard_value": entry.get("standard_value"),
                "standard_units": entry.get("standard_units"),
            }
        )
    return rows


def _normalise_string_column(
    series: pd.Series, *, uppercase: bool = False
) -> pd.Series:
    coerced = series.astype("string").str.strip()
    coerced = coerced.mask(coerced == "", other=pd.NA)
    if uppercase:
        coerced = coerced.str.upper()
    return coerced


__all__ = [
    "ActivityData",
    "ACTIVITY_COLUMN_ORDER",
    "fetch_normalize_activity",
    "generate_activity_reports",
    "normalize_activity_frame",
]

