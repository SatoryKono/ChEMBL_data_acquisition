"""Test configuration fixtures.

This module provides common pytest fixtures used across the test suite.
"""

from __future__ import annotations

import pytest

from library.config import ApiCfg, Config


@pytest.fixture(autouse=True)
def disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Disallow outbound HTTP requests during tests."""

    import requests

    def deny(*args: object, **kwargs: object) -> None:
        raise AssertionError("External network access is disabled during tests")

    monkeypatch.setattr(requests.sessions.Session, "request", deny)


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
