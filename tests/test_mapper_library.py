"""Tests for :mod:`library.mapper_library`."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, cast

import pytest
import requests

from library.config import UniprotMappingCfg
from library.mapper_library import map_chembl_to_uniprot


@dataclass
class DummyResponse:
    """Minimal stand-in for :class:`requests.Response`."""

    status_code: int
    json_data: dict[str, Any]
    reason: str = "OK"

    def json(self) -> dict[str, Any]:
        return self.json_data

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(response=self)

    def __enter__(self) -> DummyResponse:
        return self

    def __exit__(self, *exc: object) -> None:  # pragma: no cover - context cleanup
        return None


class DummySession:
    """Session returning predefined responses and recording timeouts."""

    def __init__(self, responders: list[Callable[[], DummyResponse]]) -> None:
        self._responders = responders
        self.timeouts: list[float | tuple[float, float] | None] = []

    def post(
        self, url: str, data: Any | None = None, timeout: float | None = None
    ) -> DummyResponse:
        self.timeouts.append(timeout)
        return self._responders.pop(0)()

    def get(self, url: str, timeout: float | None = None) -> DummyResponse:
        self.timeouts.append(timeout)
        return self._responders.pop(0)()


def test_map_chembl_to_uniprot_uses_config_timeout() -> None:
    """Ensure network requests use ``cfg.timeout`` instead of a hard-coded value."""

    def _run_response() -> DummyResponse:
        return DummyResponse(200, {"jobId": "1"})

    def _status_run() -> DummyResponse:
        return DummyResponse(200, {"jobStatus": "RUNNING"})

    def _status_finish() -> DummyResponse:
        return DummyResponse(200, {"jobStatus": "FINISHED"})

    def _result_response() -> DummyResponse:
        return DummyResponse(200, {"results": [{"to": {"primaryAccession": "P12345"}}]})

    session = DummySession(
        [
            _run_response,
            _status_run,
            _status_finish,
            _result_response,
        ]
    )

    cfg = UniprotMappingCfg(poll_interval=0.1, timeout=5.0)
    assert (
        map_chembl_to_uniprot("CHEMBL1", cfg, cast(requests.Session, session))
        == "P12345"
    )
    assert session.timeouts and all(t == cfg.timeout for t in session.timeouts)


def test_map_chembl_to_uniprot_handles_server_error() -> None:
    """5xx responses are reported as :class:`ValueError`."""

    def _run_error() -> DummyResponse:
        return DummyResponse(500, {})

    session = DummySession([_run_error])
    cfg = UniprotMappingCfg(poll_interval=0.1, timeout=5.0)
    with pytest.raises(ValueError):
        map_chembl_to_uniprot("CHEMBL1", cfg, cast(requests.Session, session))
