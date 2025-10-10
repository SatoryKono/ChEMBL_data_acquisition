"""Minimal ChEMBL client helpers for tissue metadata."""

from __future__ import annotations

from collections.abc import Iterable
from urllib.parse import urljoin

import pandas as pd

from ...clients import ChemblClient, _chunked
from ...common.log import logger
from ...common.pandas_utils import json_normalize_pyarrow
from ...config import ApiCfg

TISSUE_BASE_COLUMNS: list[str] = [
    "tissue_chembl_id",
    "pref_name",
    "uberon_id",
    "efo_id",
    "bto_id",
    "caloha_id",
]
"""Canonical column order for tissue records retrieved from ChEMBL."""

TISSUE_COLUMN_ORDER: list[str] = [
    *TISSUE_BASE_COLUMNS,
    "pipeline_version",
    "timestamp_utc",
]
"""Final column order after pipeline metadata enrichment."""


def _normalise_records(frame: pd.DataFrame) -> pd.DataFrame:
    """Ensure that ``frame`` exposes the expected columns in stable order."""

    work = frame.copy()
    for column in TISSUE_BASE_COLUMNS:
        if column not in work.columns:
            work[column] = pd.NA
    return work.reindex(columns=TISSUE_BASE_COLUMNS)


def get_tissues(
    ids: Iterable[str],
    *,
    cfg: ApiCfg,
    client: ChemblClient,
    chunk_size: int = 20,
    timeout: float | None = None,
) -> pd.DataFrame:
    """Fetch tissue metadata for ``ids`` from the ChEMBL API."""

    valid = [identifier for identifier in ids if identifier not in {"", "#N/A"}]
    if not valid:
        return pd.DataFrame(columns=TISSUE_BASE_COLUMNS)

    records: list[pd.DataFrame] = []
    base = f"{cfg.chembl_base.rstrip('/')}/tissue.json?format=json"
    effective_timeout = timeout if timeout is not None else cfg.timeout_read

    for chunk in _chunked(valid, max(1, chunk_size)):
        chunk_key = ",".join(chunk)
        logger.info(
            "chunk_start", extra={"stage": "chunk_start", "chunk_key": chunk_key}
        )
        url = f"{base}&tissue_chembl_id__in={chunk_key}&limit={len(chunk)}"
        chunk_frames: list[pd.DataFrame] = []
        next_url: str | None = url
        while next_url:
            data = client.request_json(next_url, cfg=cfg, timeout=effective_timeout)
            items = data.get("tissues") or data.get("tissue") or []
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
            logger.info(
                "chunk_done", extra={"stage": "chunk_done", "chunk_key": chunk_key}
            )
        else:
            logger.info(
                "chunk_skip", extra={"stage": "chunk_skip", "chunk_key": chunk_key}
            )

    if not records:
        return pd.DataFrame(columns=TISSUE_BASE_COLUMNS)

    df = pd.concat(records, ignore_index=True, sort=False)
    df = _normalise_records(df)
    return df.reset_index(drop=True)


__all__ = ["TISSUE_BASE_COLUMNS", "TISSUE_COLUMN_ORDER", "get_tissues"]
