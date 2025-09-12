"""Tests for command module wrappers in :mod:`library.cli.commands`."""

from __future__ import annotations

import importlib

import pytest

COMMANDS = [
    "get_activity_data",
    "get_assay_data",
    "get_target_data",
    "get_document_data",
    "get_document_type",
    "get_input_initialisation",
    "get_testitem_data",
    "csv_utils",
    "mapper",
    "table_quality",
    "chunk_io",
    "get_activities",
    "check_determinism",
]


@pytest.mark.parametrize("command", COMMANDS)
def test_command_import_and_run(monkeypatch: pytest.MonkeyPatch, command: str) -> None:
    """Each command module should be importable and runnable.

    The actual script invocation is patched to avoid executing heavy logic.
    """

    module = importlib.import_module(f"library.cli.commands.{command}")
    monkeypatch.setattr(module, "_run", lambda mod, argv=None: 0)
    result = module.main([])
    assert result == 0
