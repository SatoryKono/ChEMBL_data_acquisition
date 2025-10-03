"""Tests for the Python version check utility."""

from __future__ import annotations

import sys
from types import SimpleNamespace

import pytest

from library.version import require_python_version


def _mock_version(major: int, minor: int) -> object:
    """Return a minimal ``sys.version_info`` mock."""
    return SimpleNamespace(
        major=major, minor=minor, micro=0, releaselevel="final", serial=0
    )


def test_require_python_version_pass(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "version_info", _mock_version(3, 12))
    require_python_version((3, 12))


def test_require_python_version_fail(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "version_info", _mock_version(3, 11))
    with pytest.raises(RuntimeError):
        require_python_version((3, 12))
