"""Runtime helpers for configuration-dependent services."""

from __future__ import annotations

from typing import Any

from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ..common.rate_limiter import configure_limiter_cache
from .models import ApiCfg, Config, CrossRefCfg, OpenAlexCfg, RetryCfg

__all__ = [
    "configure_rate_limiters",
    "session_with_retry",
    "openalex_session",
    "crossref_session",
]


def configure_rate_limiters(cfg: Config) -> None:
    """Configure global rate limiter cache using *cfg*."""

    configure_limiter_cache(cfg.rate.limiter_cache_maxsize, cfg.rate.limiter_cache_ttl)


def session_with_retry(api: ApiCfg, retry: RetryCfg) -> Session:
    """Return an HTTP session configured for retries and user agent."""

    session = Session()
    # Disable urllib3-level retries to avoid duplicate retry loops and noisy
    # connection pool warnings.  All HTTP retry logic is implemented by the
    # higher-level clients (for example :class:`ChemblClient`) where the
    # behaviour is deterministic and centrally logged.  Leaving retries
    # enabled in ``HTTPAdapter`` leads to retries that we cannot observe or
    # control, and in offline environments this manifests as repeated
    # ``NameResolutionError`` warnings from urllib3.  Setting every retry
    # counter to ``0`` keeps adapter configuration predictable while still
    # exposing metadata such as ``backoff_max`` for introspection in tests.
    retry_kwargs: dict[str, Any] = {
        "total": 0,
        "connect": 0,
        "read": 0,
        "redirect": 0,
        "status": 0,
        "other": 0,
        "backoff_factor": retry.backoff_factor,
        "status_forcelist": retry.status_forcelist,
        "allowed_methods": None,
        "raise_on_status": False,
        "raise_on_redirect": False,
    }
    if retry.backoff_cap is not None:
        retry_kwargs["backoff_max"] = retry.backoff_cap

    retry_cfg = Retry(**retry_kwargs)
    adapter = HTTPAdapter(max_retries=retry_cfg)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers["User-Agent"] = api.user_agent
    return session


def _session_with_mailto_header(
    api: ApiCfg, retry: RetryCfg, mailto: str, verify: bool | str
) -> Session:
    session = session_with_retry(api, retry)
    session.headers["mailto"] = mailto
    session.verify = verify
    return session


def openalex_session(api: ApiCfg, retry: RetryCfg, cfg: OpenAlexCfg) -> Session:
    """Return a session configured for OpenAlex requests."""

    return _session_with_mailto_header(api, retry, cfg.mailto, cfg.verify)


def crossref_session(api: ApiCfg, retry: RetryCfg, cfg: CrossRefCfg) -> Session:
    """Return a session configured for CrossRef requests."""

    return _session_with_mailto_header(api, retry, cfg.mailto, cfg.verify)
