"""Tests for :mod:`library.utils.bootstrap`."""

from __future__ import annotations

import importlib
import sys
from types import ModuleType

import pytest

from library.utils import bootstrap


def test_ensure_project_root_supports_importing_library_schemas(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The helper pre-pends the repository root and enables ``library.schemas``."""

    project_root = str(bootstrap._PROJECT_ROOT)
    sanitized_path = [entry for entry in sys.path if entry != project_root]
    monkeypatch.setattr(sys, "path", sanitized_path, raising=False)
    monkeypatch.delitem(sys.modules, "library.schemas", raising=False)

    bootstrap.ensure_project_root()

    assert sys.path[0] == project_root

    module: ModuleType = importlib.import_module("library.schemas")
    assert module.__name__ == "library.schemas"
