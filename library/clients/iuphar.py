"""HTTP and file-loading utilities for interacting with the IUPHAR service.

The helpers in this module provide a thin client layer that encapsulates all
network and CSV loading concerns for the rest of the IUPHAR integration.  They
mirror the retry, rate limiting and column validation behaviour from the
original :mod:`library.iuphar_library` implementation while keeping the data
transformation logic decoupled from IO primitives.
"""

from __future__ import annotations

import io
import random
import threading
from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote

import pandas as pd
import requests
from requests import Session

from ..config import ApiCfg, IupharCfg, RetryCfg, session_with_retry
from ..log import logger
from ..rate_limiter import get_limiter, sleep

__all__ = [
    "EXPECTED_FAMILY_COLUMNS",
    "EXPECTED_TARGET_COLUMNS",
    "download_gtp_to_hgnc_mapping",
    "download_gtp_to_uniprot_mapping",
    "init_session",
    "load_families",
    "load_targets",
    "query_gene_symbol",
]

# Default session uses the library-wide API configuration; callers should
# override via :func:`init_session` when custom settings are required.
_session: Session = session_with_retry(ApiCfg(), RetryCfg())
_session_lock = threading.Lock()

EXPECTED_TARGET_COLUMNS: tuple[str, ...] = (
    "target_id",
    "uniprot_id",
    "hgnc_name",
    "hgnc_id",
    "gene_name",
    "synonyms",
    "family_id",
    "target_name",
)

EXPECTED_FAMILY_COLUMNS: tuple[str, ...] = (
    "family_id",
    "family_name",
    "parent_family_id",
    "target_id",
    "type",
)


def init_session(api: ApiCfg, retry: RetryCfg) -> None:
    """Initialise the shared HTTP session used by the IUPHAR client."""

    global _session
    old_session: Session | None = None
    with _session_lock:
        old_session = _session
        _session = session_with_retry(api, retry)

    if old_session is not None and old_session is not _session:
        old_session.close()


def _validate_columns(df: pd.DataFrame, expected: Iterable[str]) -> None:
    """Validate that *df* contains the *expected* columns."""

    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")


def load_targets(path: str | Path, *, encoding: str = "utf-8") -> pd.DataFrame:
    """Load the ``_IUPHAR_target.csv`` file."""

    df = pd.read_csv(path, encoding=encoding, dtype=str).fillna("")
    # Normalise legacy column names to the expected lowercase form
    df = df.rename(
        columns={
            "GuidetoPHARMACOLOGY": "target_id",
            "HGNC_NAME": "hgnc_name",
            "HGNC_name": "hgnc_name",
            "HGNC_ID": "hgnc_id",
            "HGNC_id": "hgnc_id",
            "swissprot": "uniprot_id",
        }
    )
    _validate_columns(df, EXPECTED_TARGET_COLUMNS)
    return df


def load_families(path: str | Path, *, encoding: str = "utf-8") -> pd.DataFrame:
    """Load the ``_IUPHAR_family.csv`` file."""

    df = pd.read_csv(path, encoding=encoding, dtype=str).fillna("")
    _validate_columns(df, EXPECTED_FAMILY_COLUMNS)
    return df


def _download_csv(url: str, cfg: IupharCfg, retry: RetryCfg) -> pd.DataFrame:
    """Download a CSV file from *url* with retry handling."""

    timeout = (cfg.timeout_connect, cfg.timeout_read)
    limiter = get_limiter("iuphar", cfg.rps, cfg.burst)

    for attempt in range(1, retry.max_attempts + 1):
        limiter.acquire()
        try:
            with _session_lock:
                with _session.get(url, timeout=timeout) as resp:
                    resp.raise_for_status()
                    return pd.read_csv(io.StringIO(resp.text))
        except requests.RequestException as exc:  # pragma: no cover - network
            if attempt >= retry.max_attempts:
                logger.error(
                    "request_error",
                    url=url,
                    attempt=attempt,
                    error=str(exc),
                )
                raise
            backoff = retry.backoff_factor * (2 ** (attempt - 1))
            jitter = random.uniform(0, backoff)
            sleep(backoff + jitter)

    raise RuntimeError("Failed to download CSV")


def _data_base(cfg: IupharCfg) -> str:
    """Return the DATA base URL derived from *cfg*."""

    base_root = cfg.base.rstrip("/")
    if base_root.endswith("services"):
        base_root = base_root.rsplit("/", 1)[0]
    return f"{base_root}/DATA"


def download_gtp_to_uniprot_mapping(cfg: IupharCfg, retry: RetryCfg) -> pd.DataFrame:
    """Download the GtP to UniProt mapping CSV."""

    data_base = _data_base(cfg)
    url = f"{data_base}/GtP_to_UniProt_mapping.csv"
    return _download_csv(url, cfg, retry)


def download_gtp_to_hgnc_mapping(cfg: IupharCfg, retry: RetryCfg) -> pd.DataFrame:
    """Download the GtP to HGNC mapping CSV."""

    data_base = _data_base(cfg)
    url = f"{data_base}/GtP_to_HGNC_mapping.csv"
    return _download_csv(url, cfg, retry)


def _query_gene_symbol(
    gene_name: str, cfg: IupharCfg, retry: RetryCfg
) -> dict[str, Any]:
    """Return the first IUPHAR result for *gene_name*."""

    base = cfg.base.rstrip("/")
    url = f"{base}/targets/?geneSymbol={quote(gene_name)}"
    timeout = (cfg.timeout_connect, cfg.timeout_read)
    limiter = get_limiter("iuphar", cfg.rps, cfg.burst)

    for attempt in range(1, retry.max_attempts + 1):
        limiter.acquire()
        try:
            with _session_lock:
                with _session.get(url, timeout=timeout) as response:
                    response.raise_for_status()
                    data = cast(list[dict[str, Any]], response.json())
                    return data[0] if data else {}
        except requests.RequestException as exc:  # pragma: no cover - network errors
            if attempt >= retry.max_attempts:
                logger.error(
                    "request_error",
                    url=url,
                    attempt=attempt,
                    error=str(exc),
                )
                break
            backoff = retry.backoff_factor * (2 ** (attempt - 1))
            jitter = random.uniform(0, backoff)
            sleep(backoff + jitter)
    return {}


def query_gene_symbol(
    gene_name: str, cfg: IupharCfg, retry: RetryCfg
) -> dict[str, Any]:
    """Public wrapper around :func:`_query_gene_symbol`."""

    return _query_gene_symbol(gene_name, cfg, retry)
