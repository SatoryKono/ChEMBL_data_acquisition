"""Shared HTTP utilities for ChEMBL API access."""

from __future__ import annotations

from typing import Any, Iterable, Iterator

import logging

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

logger = logging.getLogger(__name__)

# Configure session with retry/backoff for robustness
_retry = Retry(
    total=3,
    backoff_factor=1.0,
    status_forcelist=[500, 502, 503, 504],
    allowed_methods=["GET"],
)
_session = requests.Session()
_adapter = HTTPAdapter(max_retries=_retry)
_session.mount("http://", _adapter)
_session.mount("https://", _adapter)


def request_json(url: str, *, timeout: float = 30.0) -> dict[str, Any]:
    """Return JSON content from *url*.

    Parameters
    ----------
    url:
        API endpoint to query.
    timeout:
        Maximum number of seconds to wait for the response.

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
    with _session.get(url, timeout=timeout) as response:
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
