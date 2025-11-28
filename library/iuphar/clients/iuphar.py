"""HTTP and file-loading utilities for interacting with the IUPHAR service.

The helpers in this module provide a thin client layer that encapsulates all
network and CSV loading concerns for the rest of the IUPHAR integration.  They
mirror the retry, rate limiting and column validation behaviour from the
original :mod:`library.iuphar.integration.iuphar_library` implementation while keeping the data
transformation logic decoupled from IO primitives.
"""

from __future__ import annotations

import io
import random
import threading
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, cast
from urllib.parse import quote

import pandas as pd
import requests
from requests import Session

from ..common.log import logger
from ..common.rate_limiter import get_limiter, sleep
from ..config.models import ApiCfg, IupharCfg, RetryCfg
from ..config.runtime import session_with_retry

__all__ = [
    "EXPECTED_FAMILY_COLUMNS",
    "EXPECTED_TARGET_COLUMNS",
    "download_gtp_to_hgnc_mapping",
    "download_gtp_to_uniprot_mapping",
    "init_session",
    "get_session",
    "load_families",
    "load_targets",
    "query_gene_symbol",
]

# Default session uses the library-wide API configuration; callers should
# override via :func:`init_session` when custom settings are required.
_session_lock = threading.RLock()
_session_condition = threading.Condition(_session_lock)


def _make_session_factory(api: ApiCfg, retry: RetryCfg) -> Callable[[], Session]:
    def factory() -> Session:
        return session_with_retry(api, retry)

    return factory


_session_factory: Callable[[], Session] = _make_session_factory(ApiCfg(), RetryCfg())
_session_local = threading.local()
_sessions: set[Session] = set()
_active_requests = 0
_jitter_provider: Callable[[float], float] | None = None


def _close_sessions(sessions: Iterable[Session]) -> None:
    for session in sessions:
        close = getattr(session, "close", None)
        if callable(close):
            close()


def _get_or_create_session_locked() -> Session:
    session_attr = getattr(_session_local, "session", None)
    session = session_attr if isinstance(session_attr, Session) else None
    if session is None:
        session = _session_factory()
        _sessions.add(session)
        _session_local.session = session
    return session


@contextmanager
def _session_context() -> Iterator[Session]:
    global _active_requests
    with _session_condition:
        session = _get_or_create_session_locked()
        _active_requests += 1
    try:
        yield session
    finally:
        with _session_condition:
            _active_requests -= 1
            if _active_requests == 0:
                _session_condition.notify_all()


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


def init_session(
    api: ApiCfg,
    retry: RetryCfg,
    *,
    jitter: Callable[[float], float] | None = None,
) -> None:
    """Initialise the shared HTTP session used by the IUPHAR client."""

    global _session_factory, _session_local, _jitter_provider

    factory = _make_session_factory(api, retry)
    with _session_condition:
        while _active_requests > 0:
            _session_condition.wait()
        old_sessions = list(_sessions)
        _sessions.clear()
        _session_factory = factory
        _session_local = threading.local()
        _jitter_provider = jitter if jitter is not None else retry.build_jitter()

    _close_sessions(old_sessions)


def get_session() -> Session:
    """Expose the configured :class:`requests.Session` instance."""

    with _session_condition:
        return _get_or_create_session_locked()


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

    with _session_condition:
        jitter_provider = _jitter_provider

    for attempt in range(1, retry.max_attempts + 1):
        limiter.acquire()
        try:
            with _session_context() as session:
                with session.get(url, timeout=timeout) as resp:
                    resp.raise_for_status()
                    csv_buffer = io.StringIO(resp.text)
                    return pd.read_csv(csv_buffer, dtype=str).fillna("")
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
            if backoff > 0:
                if jitter_provider is not None:
                    jitter_value = jitter_provider(backoff)
                else:
                    jitter_value = random.uniform(0, backoff)
            else:
                jitter_value = 0.0
            sleep(backoff + jitter_value)

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

    with _session_condition:
        jitter_provider = _jitter_provider

    for attempt in range(1, retry.max_attempts + 1):
        limiter.acquire()
        try:
            with _session_context() as session:
                with session.get(url, timeout=timeout) as response:
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
            if backoff > 0:
                if jitter_provider is not None:
                    jitter_value = jitter_provider(backoff)
                else:
                    jitter_value = random.uniform(0, backoff)
            else:
                jitter_value = 0.0
            sleep(backoff + jitter_value)
    return {}


def query_gene_symbol(
    gene_name: str, cfg: IupharCfg, retry: RetryCfg
) -> dict[str, Any]:
    """Public wrapper around :func:`_query_gene_symbol`."""

    return _query_gene_symbol(gene_name, cfg, retry)
