"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

from typing import Any, Iterable, Iterator

import logging

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import Config

logger = logging.getLogger(__name__)


def request_json(
    cfg: Config, url: str, *, timeout: float | None = None
) -> dict[str, Any]:
    """Return JSON content from ``url``.

    Parameters
    ----------
    cfg:
        Application configuration containing timeout and retry settings.
    url:
        API endpoint to query.
    timeout:
        Maximum number of seconds to wait for the response. When ``None`` the
        value from :attr:`cfg.timeouts.read` is used.

    Returns
    -------
    dict[str, Any]
        Parsed JSON document.

    Raises
    ------
    requests.RequestException
        If the HTTP request fails.
    ValueError
        If the response body is not valid JSON.

    Notes
    -----
    A new :class:`requests.Session` is created for each call with retry
    behaviour configured according to ``cfg``. This keeps the function
    stateless and avoids hidden global configuration.
    """

    retry = Retry(
        total=cfg.rate_limits.max_retries,
        backoff_factor=cfg.rate_limits.backoff_factor,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    session = requests.Session()
    session.mount("http://", adapter)
    session.mount("https://", adapter)

    effective_timeout = timeout if timeout is not None else cfg.timeouts.read
    with session.get(url, timeout=effective_timeout) as response:
        response.raise_for_status()
        return response.json()


def _chunked(items: Iterable[str], size: int) -> Iterator[list[str]]:
    """Yield ``size``-sized lists from *items*.

    Parameters
    ----------
    items:
        Iterable of identifiers to split.
    size:
        Desired chunk size; must be positive.

    Yields
    ------
    list[str]
        Subsequences of ``items`` with at most ``size`` elements.

    Raises
    ------
    ValueError
        If ``size`` is not a positive integer.

    """
    if size <= 0:
        raise ValueError("size must be a positive integer")

    chunk: list[str] = []
    for item in items:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk
