"""ChEMBL client helpers for retrieving cell line metadata."""

from __future__ import annotations

from collections.abc import Iterable
from urllib.parse import urljoin

import pandas as pd

from ...clients import ChemblClient, _chunked
from ...common.log import logger
from ...common.pandas_utils import json_normalize_pyarrow
from ...config import ApiCfg

MAX_CELLLINE_CHUNK_SIZE = 1000

CELL_LINE_BASE_COLUMNS: list[str] = [
    "cell_chembl_id",
    "cell_name",
    "cell_description",
    "cell_id",
    "cell_source_organism",
    "cell_source_tax_id",
    "cell_source_tissue",
    "cellosaurus_id",
    "cl_lincs_id",
    "clo_id",
    "efo_id",
]
"""Canonical column order for cell line records retrieved from ChEMBL."""

CELL_LINE_COLUMN_ORDER: list[str] = [
    *CELL_LINE_BASE_COLUMNS,
    "pipeline_version",
    "timestamp_utc",
]
"""Final column order after pipeline metadata enrichment."""


def _ensure_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` with all :data:`CELL_LINE_BASE_COLUMNS` present."""

    result = frame.copy()
    for column in CELL_LINE_BASE_COLUMNS:
        if column not in result.columns:
            result[column] = pd.NA
    return result.reindex(columns=CELL_LINE_BASE_COLUMNS)


def _coerce_types(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of ``frame`` with stable nullable dtypes."""

    result = frame.copy()
    string_columns = [
        "cell_chembl_id",
        "cell_name",
        "cell_description",
        "cell_source_organism",
        "cell_source_tissue",
        "cellosaurus_id",
        "cl_lincs_id",
        "clo_id",
        "efo_id",
    ]
    integer_columns = ["cell_id", "cell_source_tax_id"]

    for column in string_columns:
        result[column] = result[column].astype(pd.StringDtype())

    for column in integer_columns:
        if column in result.columns:
            result[column] = result[column].astype("Int64")
        else:
            result[column] = pd.Series(pd.NA, dtype="Int64")

    return result


def _normalise_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure cell line frames expose the expected layout and dtypes."""

    work = _ensure_columns(frame)
    work = work.replace({"": pd.NA})
    work = _coerce_types(work)
    return work


def get_cell_lines(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 20,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch cell line metadata for ``ids`` from the ChEMBL API."""

    valid = [identifier for identifier in ids if identifier not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=CELL_LINE_BASE_COLUMNS)

    base = f"{cfg.chembl_base.rstrip('/')}/cell_line.json?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read

    records: list[pd.DataFrame] = []
    effective_chunk_size = max(1, min(int(chunk_size), MAX_CELLLINE_CHUNK_SIZE))
    if effective_chunk_size != chunk_size:
        logger.debug(
            "cellline_chunk_clamped",
            extra={
                "requested_chunk_size": chunk_size,
                "effective_chunk_size": effective_chunk_size,
                "limit": MAX_CELLLINE_CHUNK_SIZE,
                "stage": "chunk_prepare",
            },
        )

    for chunk in _chunked(valid, effective_chunk_size):
        chunk_key = ",".join(chunk)
        logger.info("cellline_chunk_start", chunk_key=chunk_key)
        url = f"{base}&cell_chembl_id__in={chunk_key}&limit={len(chunk)}"
        chunk_frames: list[pd.DataFrame] = []
        next_url: str | None = url
        while next_url:
            data = client.request_json(next_url, cfg=cfg, timeout=effective_timeout)
            items = data.get("cell_lines") or data.get("cell_line") or []
            if items:
                df_chunk = json_normalize_pyarrow(items).dropna(
                    axis="columns", how="all"
                )
                if not df_chunk.empty:
                    chunk_frames.append(_normalise_records(df_chunk))
            page_meta = data.get("page_meta") or {}
            next_token = page_meta.get("next")
            next_url = urljoin(cfg.chembl_base, next_token) if next_token else None
        if chunk_frames:
            records.append(pd.concat(chunk_frames, ignore_index=True, sort=False))
            logger.info("cellline_chunk_done", chunk_key=chunk_key)
        else:
            logger.info("cellline_chunk_skip", chunk_key=chunk_key)

    if not records:
        return pd.DataFrame(columns=CELL_LINE_BASE_COLUMNS)

    df = pd.concat(records, ignore_index=True, sort=False)
    df = _normalise_records(df)
    return df.reset_index(drop=True)


__all__ = [
    "CELL_LINE_BASE_COLUMNS",
    "CELL_LINE_COLUMN_ORDER",
    "get_cell_lines",
]
