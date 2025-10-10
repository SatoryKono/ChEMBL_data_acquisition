"""Smoke tests ensuring CLI compatibility wrappers stay importable."""

from __future__ import annotations

import importlib

import pytest


@pytest.mark.unit
def test_get_document_data_cli__import_smoke() -> None:
    """Compatibility wrapper must remain importable and proxy the CLI module."""

    cli_module = importlib.import_module("scripts.get_document_data")
    impl_module = importlib.import_module("library.cli.commands.get_document_data")

    assert cli_module is impl_module
    assert hasattr(cli_module, "main")
