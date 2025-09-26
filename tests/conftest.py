"""Test configuration fixtures.

This module provides common pytest fixtures used across the test suite.
"""

from __future__ import annotations

import pytest

from library.config import Config


@pytest.fixture(autouse=True)
def _disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Prevent external HTTP requests during tests for determinism."""

    try:
        import requests
    except ModuleNotFoundError:  # pragma: no cover - requests optional in env
        return

    def _deny_request(self, method, url, *args, **kwargs):  # type: ignore[override]
        msg = (
            "External HTTP requests are disabled during tests; "
            f"attempted {method} {url}"
        )
        raise RuntimeError(msg)

    monkeypatch.setattr("requests.sessions.Session.request", _deny_request)


@pytest.fixture()
def cfg() -> Config:
    """Return a baseline :class:`~library.config.Config` instance for tests.

    The configuration requires a valid ``api.user_agent`` value; the fixture
    supplies a deterministic placeholder address so that individual tests do
    not need to construct :class:`Config` instances manually.
    """

    cfg = Config()
    cfg.api.user_agent = "test@example.org"
    return cfg


@pytest.fixture()
def duplicate_document_ids() -> list[str]:
    """Return sample document IDs including duplicates for testing."""

    return ["CHEMBL1", "CHEMBL1", "CHEMBL2"]
