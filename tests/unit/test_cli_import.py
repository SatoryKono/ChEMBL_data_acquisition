"""Smoke tests ensuring CLI compatibility wrappers stay importable."""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.unit
def test_get_document_data_cli__import_smoke() -> None:
    """Compatibility wrapper must remain importable and proxy the CLI module."""

    cli_module = importlib.import_module("scripts.get_document_data")
    impl_module = importlib.import_module("library.cli.commands.get_document_data")

    assert cli_module is impl_module
    assert hasattr(cli_module, "main")


@pytest.mark.unit
def test_get_document_data_cli__script_help() -> None:
    """Invoking the legacy script entry point must surface the CLI help."""

    script_path = Path("scripts/get_document_data.py")
    result = subprocess.run(
        [sys.executable, str(script_path), "--help"],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0
    combined_output = f"{result.stdout}\n{result.stderr}".lower()
    assert "usage" in combined_output
