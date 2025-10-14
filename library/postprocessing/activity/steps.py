"""Postprocessing helpers for deterministic ChEMBL activity exports."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

from library.clients.chembl_client import ChemblClient
from library.schemas.activity_schema import ACTIVITY_COLUMN_ORDER, validate_activity_records
from library.utils.logging import Logger, get_logger
from library.utils.qc_report import TableQualityProfiler, build_reports_from_profiler

DEFAULT_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
_DEFAULT_PAGE_SIZE = 500
_ACTIVITY_FIELDS: tuple[str, ...] = ACTIVITY_COLUMN_ORDER
_STRING_COLUMNS: tuple[str, ...] = (
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
    "standard_units",
)
_UPPERCASE_COLUMNS: tuple[str, ...] = (
    "assay_chembl_id",
    "target_chembl_id",
    "standard_relation",
    "standard_units",
    "standard_type",
)
_NUMERIC_COLUMNS: tuple[str, ...] = ("standard_value",)
ACTIVITY_SORT_COLUMNS: tuple[str, ...] = (
    "activity_id",
    "assay_chembl_id",
    "target_chembl_id",
    "standard_type",
    "standard_relation",
)

__all__ = [
    "ACTIVITY_COLUMN_ORDER",
    "ACTIVITY_SORT_COLUMNS",
    "build_activity_reports",
    "fetch_normalize_activity",
    "normalize_activity_frame",
    "run_activity_pipeline",
]


def _prepare_client(client: ChemblClient | None, *, base_url: str | None) -> ChemblClient:
    """Return a configured :class:`ChemblClient`."""

    if client is not None:
        return client
    resolved_base = base_url or DEFAULT_BASE_URL
    return ChemblClient(base_url=resolved_base)


def _fetch_page(
    client: ChemblClient,
    *,
    params: Mapping[str, Any],
) -> Mapping[str, Any]:
    """Return the JSON payload for the activities endpoint."""

    kwargs: dict[str, Any] = {"params": params}
    timeout = getattr(client, "timeout", None)
    if timeout is not None:
        kwargs["timeout"] = timeout
    try:
        return client.list_activities(**kwargs)
    except TypeError:  # pragma: no cover - non-standard client signature
        kwargs.pop("timeout", None)
        return client.list_activities(**kwargs)  # type: ignore[misc]


def _coerce_string(series: pd.Series) -> pd.Series:
    """Return ``series`` coerced to ``StringDtype`` with trimmed values."""

    coerced = series.astype("string").str.strip()
    return coerced.replace({"": pd.NA})


def _uppercase(series: pd.Series) -> pd.Series:
    """Return ``series`` with upper-cased textual content."""

    return series.str.upper()


def normalize_activity_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Apply deterministic coercions and schema validation to ``df``."""

    prepared = df.copy(deep=True)
    for column in _ACTIVITY_FIELDS:
        if column not in prepared.columns:
            prepared[column] = pd.NA

    for column in _STRING_COLUMNS:
        prepared[column] = _coerce_string(prepared[column])
        if column in _UPPERCASE_COLUMNS:
            prepared[column] = _uppercase(prepared[column])

    activity_numeric = pd.to_numeric(prepared["activity_id"], errors="coerce")
    prepared["activity_id"] = pd.Series(
        pd.array(activity_numeric, dtype="Int64"), index=prepared.index
    )

    for column in _NUMERIC_COLUMNS:
        numeric = pd.to_numeric(prepared[column], errors="coerce")
        prepared[column] = pd.Series(
            pd.array(numeric, dtype="Float64"), index=prepared.index
        )

    ordered = prepared.loc[:, ACTIVITY_COLUMN_ORDER]
    ordered = ordered.sort_values(
        list(ACTIVITY_SORT_COLUMNS), kind="mergesort", na_position="last"
    )
    validated = validate_activity_records(ordered.reset_index(drop=True))
    return validated


def fetch_normalize_activity(
    limit: int | None = None,
    *,
    client: ChemblClient | None = None,
    base_url: str | None = None,
    page_size: int = _DEFAULT_PAGE_SIZE,
    logger: Logger | None = None,
) -> pd.DataFrame:
    """Fetch activities from ChEMBL and return a normalized DataFrame."""

    if limit is not None and limit < 0:
        raise ValueError("limit must be non-negative")

    effective_logger = logger or get_logger(__name__).bind(stage="activity_fetch")
    effective_client = _prepare_client(client, base_url=base_url)
    owns_client = client is None

    records: list[dict[str, Any]] = []
    fetched = 0
    offset = 0
    try:
        while True:
            if limit is not None and fetched >= limit:
                break
            batch_size = page_size
            if limit is not None:
                remaining = limit - fetched
                if remaining <= 0:
                    break
                batch_size = min(batch_size, remaining)
            params = {
                "limit": batch_size,
                "offset": offset,
                "format": "json",
                "order_by": "activity_id",
                "fields": ",".join(_ACTIVITY_FIELDS),
            }
            effective_logger.info("activity_fetch_page_start", offset=offset, limit=batch_size)
            payload = _fetch_page(effective_client, params=params)
            activities = payload.get("activities", [])
            if not isinstance(activities, list):
                activities = []
            effective_logger.info(
                "activity_fetch_page_complete", offset=offset, rows=len(activities)
            )
            for entry in activities:
                if not isinstance(entry, Mapping):
                    continue
                record = {field: entry.get(field) for field in _ACTIVITY_FIELDS}
                records.append(record)
                fetched += 1
                if limit is not None and fetched >= limit:
                    break
            if not activities:
                break
            meta = payload.get("page_meta")
            has_next = False
            if isinstance(meta, Mapping):
                has_next = bool(meta.get("next"))
            if not has_next:
                break
            offset += batch_size
    finally:
        if owns_client:
            close = getattr(effective_client, "close", None)
            if callable(close):
                close()

    frame = pd.DataFrame(records, columns=_ACTIVITY_FIELDS)
    effective_logger.info("activity_fetch_complete", rows=len(frame))
    return normalize_activity_frame(frame)


def build_activity_reports(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return quality and correlation reports for ``frame``."""

    if frame.empty:
        return pd.DataFrame(), pd.DataFrame()

    profiler = TableQualityProfiler()
    profiler.consume(frame)
    quality, correlation = build_reports_from_profiler(profiler)
    if not quality.empty and "column" in quality.columns:
        quality = quality.sort_values("column", kind="mergesort").reset_index(drop=True)
    if not correlation.empty:
        correlation = correlation.sort_index(axis=0).sort_index(axis=1)
    return quality, correlation


def run_activity_pipeline(
    *,
    limit: int | None = None,
    client: ChemblClient | None = None,
    base_url: str | None = None,
    page_size: int = _DEFAULT_PAGE_SIZE,
    logger: Logger | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return the normalized dataset and QC artefacts for the activity table."""

    dataset = fetch_normalize_activity(
        limit,
        client=client,
        base_url=base_url,
        page_size=page_size,
        logger=logger,
    )
    quality, correlation = build_activity_reports(dataset)
    return dataset, correlation, quality
