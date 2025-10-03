from __future__ import annotations

import importlib


def test_mapper_batch_library_imports() -> None:
    """Ensure mapper_batch_library imports without errors."""

    module = importlib.import_module("library.integration.mapper_batch_library")
    assert module.__name__ == "library.integration.mapper_batch_library"
