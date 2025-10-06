"""Runtime helpers for configuration consumers."""

from __future__ import annotations

from typing import Any

from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ..common.rate_limiter import configure_limiter_cache
from .models import ApiCfg, Config, CrossRefCfg, OpenAlexCfg, RateCfg, RetryCfg

__all__ = [
    "configure_rate_limits",
    "ensure_dirs",
    "openalex_session",
    "crossref_session",
    "session_with_retry",
]


def session_with_retry(api: ApiCfg, retry: RetryCfg) -> Session:
    session = Session()
    retry_kwargs: dict[str, Any] = {
        "total": max(0, retry.max_attempts - 1),
        "backoff_factor": retry.backoff_factor,
        "status_forcelist": retry.status_forcelist,
        "allowed_methods": None,
    }
    if retry.backoff_cap is not None:
        retry_kwargs["backoff_max"] = retry.backoff_cap
    jitter = retry.build_jitter()
    if jitter is not None:
        retry_kwargs["backoff_jitter"] = jitter
    adapter = HTTPAdapter(max_retries=Retry(**retry_kwargs))
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers["User-Agent"] = api.user_agent
    return session


def _session_with_mailto_header(api: ApiCfg, retry: RetryCfg, mailto: str) -> Session:
    session = session_with_retry(api, retry)
    session.headers["mailto"] = mailto
    return session


def openalex_session(api: ApiCfg, retry: RetryCfg, cfg: OpenAlexCfg) -> Session:
    return _session_with_mailto_header(api, retry, cfg.mailto)


def crossref_session(api: ApiCfg, retry: RetryCfg, cfg: CrossRefCfg) -> Session:
    return _session_with_mailto_header(api, retry, cfg.mailto)


def configure_rate_limits(rate: RateCfg) -> None:
    configure_limiter_cache(rate.limiter_cache_maxsize, rate.limiter_cache_ttl)


def ensure_dirs(cfg: Config) -> None:
    for path in (cfg.io.output_dir, cfg.io.cache_dir):
        if path.exists():
            if not path.is_dir():
                raise NotADirectoryError(f"{path} is not a directory")
        else:
            if cfg.io.exist_ok:
                path.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(f"{path} does not exist")
