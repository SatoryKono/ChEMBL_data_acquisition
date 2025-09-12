"""Tests for :mod:`library.mapper_library`."""

from __future__ import annotations

import io
import json
import urllib.request
from typing import Any

import pytest

import library.mapper_library as mapper_library
from library.config import UniprotMappingCfg
from library.mapper_library import map_chembl_to_uniprot


class _BytesIOContext:
    """Context manager returning a :class:`io.BytesIO` object."""

    def __init__(self, data: bytes) -> None:
        self._bio = io.BytesIO(data)

    def __enter__(self) -> io.BytesIO:
        return self._bio

    def __exit__(self, *args: object) -> None:
        self._bio.close()


def test_map_chembl_to_uniprot_uses_config_timeout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure network requests use ``cfg.timeout`` instead of a hard-coded value."""
    mapper_library._map_chembl_to_uniprot_cached.cache_clear()

    timeouts: list[float | None] = []

    def fake_urlopen(
        url: str, data: bytes | None = None, timeout: float | None = None
    ) -> _BytesIOContext:
        timeouts.append(timeout)
        if url.endswith("/run"):
            content: dict[str, Any] = {"jobId": "1"}
        elif "/status/" in url:
            if fake_urlopen.calls == 0:
                fake_urlopen.calls += 1
                content = {"jobStatus": "RUNNING"}
            else:
                content = {"jobStatus": "FINISHED"}
        else:
            content = {"results": [{"to": {"primaryAccession": "P12345"}}]}
        return _BytesIOContext(json.dumps(content).encode())

    fake_urlopen.calls = 0  # type: ignore[attr-defined]

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

    cfg = UniprotMappingCfg(poll_interval=0.1, timeout=5.0)
    assert map_chembl_to_uniprot("CHEMBL1", cfg) == "P12345"
    assert timeouts and all(t == cfg.timeout for t in timeouts)


def test_map_chembl_to_uniprot_cached(monkeypatch: pytest.MonkeyPatch) -> None:
    """Repeated IDs should trigger only one HTTP call."""

    mapper_library._map_chembl_to_uniprot_cached.cache_clear()

    calls = 0

    def fake_urlopen(
        url: str, data: bytes | None = None, timeout: float | None = None
    ) -> _BytesIOContext:
        nonlocal calls
        calls += 1
        if url.endswith("/run"):
            content: dict[str, Any] = {"jobId": "1"}
        elif "/status/" in url:
            content = {"jobStatus": "FINISHED"}
        else:
            content = {"results": [{"to": {"primaryAccession": "P12345"}}]}
        return _BytesIOContext(json.dumps(content).encode())

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)

    cfg = UniprotMappingCfg(poll_interval=0.1, timeout=5.0)
    assert map_chembl_to_uniprot("CHEMBL2", cfg) == "P12345"
    assert map_chembl_to_uniprot("CHEMBL2", cfg) == "P12345"
    # 3 calls: run, status, results
    assert calls == 3
