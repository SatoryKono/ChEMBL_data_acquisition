"""Tests for thin wrappers in :mod:`library.cli.entrypoints`."""

from __future__ import annotations

import importlib

import pytest


@pytest.mark.unit
@pytest.mark.parametrize(
    ("module_name", "expected_target"),
    [
        ("get_data", "get_data"),
        ("assay", "get_assay_data"),
        ("target", "get_target_data"),
        ("document", "get_document_data"),
        ("document_type", "get_document_type"),
        ("input_initialisation", "get_input_initialisation"),
        ("testitem", "get_testitem_data"),
        ("csv_utils", "csv_utils"),
        ("mapper", "mapper"),
        ("table_quality", "table_quality"),
        ("chunk_io", "chunk_io"),
        ("activities", "get_activities"),
        ("check_determinism", "check_determinism"),
    ],
)
def test_entrypoint_wrappers__delegate_to_legacy(
    module_name: str,
    expected_target: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure thin wrappers forward calls to :mod:`library.cli.commands`."""

    observed: dict[str, object] = {}

    def _fake_execute(target: str, argv: list[str] | None) -> int:
        observed["target"] = target
        observed["argv"] = argv
        return 0

    monkeypatch.setattr(
        "library.cli.entrypoints._legacy.execute",
        _fake_execute,
        raising=True,
    )

    module = importlib.import_module(f"library.cli.entrypoints.{module_name}")
    argv = ["--flag", "value"]

    main_func = getattr(module, "main")
    assert callable(main_func)

    exit_code = main_func(argv)

    assert exit_code == 0
    assert observed == {"target": expected_target, "argv": argv}

