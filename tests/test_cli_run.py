"""Regression tests for CLI invocation helpers."""

from __future__ import annotations

from types import SimpleNamespace

from library.cli.commands import _run


def test_run_invokes_module_main(monkeypatch) -> None:
    """Ensure ``_run`` forwards argv to the resolved module."""

    captured: dict[str, object] = {}

    def fake_import(name: str):
        captured["module"] = name
        return SimpleNamespace(main=lambda argv=None: (captured.update(argv=argv) or 5))

    monkeypatch.setattr("library.cli.commands.import_module", fake_import)

    exit_code = _run("scripts.get_activity_data", ["--print-config"])

    assert exit_code == 5
    assert captured["module"] == "scripts.get_activity_data"
    assert captured["argv"] == ["--print-config"]
