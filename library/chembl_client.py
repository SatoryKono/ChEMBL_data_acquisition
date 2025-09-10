"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

from typing import Any, Iterable, Iterator

import logging

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import ApiCfg

logger = logging.getLogger(__name__)


def request_json(
    url: str, *, cfg: ApiCfg, timeout: float | None = None
) -> dict[str, Any]:
    """Return JSON content from *url*.

    Parameters
    ----------
    url:
        API endpoint to query.
    cfg:
        Configuration providing timeout and retry settings.
    timeout:
        Optional override for the read timeout in seconds.

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

    """
    retry = Retry(
        total=cfg.retries,
        backoff_factor=cfg.backoff_factor,
        status_forcelist=[500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry)
    with requests.Session() as session:
        session.mount("http://", adapter)
        session.mount("https://", adapter)
        read_timeout = timeout if timeout is not None else cfg.timeout_read
        with session.get(url, timeout=(cfg.timeout_connect, read_timeout)) as response:
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
