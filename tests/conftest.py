"""Test configuration fixtures.

This module provides common pytest fixtures used across the test suite.
"""

from __future__ import annotations

import pytest

from library.config import ApiCfg, Config


@pytest.fixture()
def cfg() -> Config:
    """Return a baseline :class:`~library.config.Config` instance for tests.

    The configuration requires a valid ``api.user_agent`` value; the fixture
    supplies a deterministic placeholder address so that individual tests do
    not need to construct :class:`Config` instances manually.
    """

    return Config(api=ApiCfg(user_agent="test@example.org"))


@pytest.fixture()
def duplicate_document_ids() -> list[str]:
    """Return sample document IDs including duplicates for testing."""

    return ["CHEMBL1", "CHEMBL1", "CHEMBL2"]
