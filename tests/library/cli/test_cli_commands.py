"""Tests for command module wrappers in :mod:`library.cli.commands`."""

from __future__ import annotations

import importlib
from types import SimpleNamespace

COMMAND_MODULE = importlib.import_module("library.cli.commands")

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


def test_run_passes_argv_to_parameterised_main(monkeypatch: pytest.MonkeyPatch) -> None:
    """``_run`` should forward the provided arguments to ``main``."""

    captured: dict[str, object] = {}

    def fake_main(argv):  # type: ignore[no-untyped-def]
        captured["argv"] = argv
        return 42

    monkeypatch.setattr(
        COMMAND_MODULE,
        "_resolve_module",
        lambda name: SimpleNamespace(main=fake_main),
    )

    exit_code = COMMAND_MODULE._run("dummy", ["--flag"])

    assert exit_code == 42
    assert captured["argv"] == ["--flag"]


def test_run_skips_argv_for_main_without_parameters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """``_run`` should avoid passing arguments when ``main`` expects none."""

    called: dict[str, bool] = {"flag": False}

    def fake_main():  # type: ignore[no-untyped-def]
        called["flag"] = True
        return 7

    monkeypatch.setattr(
        COMMAND_MODULE,
        "_resolve_module",
        lambda name: SimpleNamespace(main=fake_main),
    )

    exit_code = COMMAND_MODULE._run("dummy", ["--ignored"])

    assert exit_code == 7
    assert called["flag"] is True
