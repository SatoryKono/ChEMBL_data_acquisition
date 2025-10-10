"""Smoke tests ensuring CLI compatibility wrappers stay importable."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest


@pytest.mark.unit
def test_get_document_data_cli__import_smoke() -> None:
    """Compatibility wrapper must remain importable and proxy the CLI module."""

    cli_module = importlib.import_module("scripts.get_document_data")
    impl_module = importlib.import_module("library.cli.commands.get_document_data")

    assert cli_module is impl_module
    assert hasattr(cli_module, "main")
    assert cli_module.DEFAULT_INPUT_NAME == impl_module.DEFAULT_INPUT_NAME
    assert cli_module._EXPORT_COLUMNS is impl_module._EXPORT_COLUMNS

    parser, _log_cfg = cli_module.build_parser()
    args = parser.parse_args(["--mode", "pubmed"])

    assert args.mode == "pubmed"
    assert args.input_csv == Path(cli_module.DEFAULT_INPUT_NAME)
