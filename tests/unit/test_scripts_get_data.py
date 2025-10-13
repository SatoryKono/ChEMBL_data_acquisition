from __future__ import annotations

import importlib
import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest


@pytest.mark.unit
def test_load_module__merge_conflict_hint(monkeypatch):
    """The compatibility wrapper should surface unresolved merge conflicts."""

    script_path = Path("scripts/get_data.py").resolve()
    module_name = "tests.scripts_get_data_conflict"
    loader = SourceFileLoader(module_name, str(script_path))
    spec = importlib.util.spec_from_loader(module_name, loader)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)

    original_import_module = importlib.import_module

    def _patched_import_module(name: str, package: str | None = None):  # type: ignore[override]
        if name == "library.cli.commands.get_data":
            raise SyntaxError(
                "invalid decimal literal",
                (
                    "library/cli_utils.py",
                    201,
                    0,
                    ">>>>>>> 9f4ae5808230279e785241afb817c7d4f744ff3a\n",
                ),
            )
        return original_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", _patched_import_module)

    try:
        with pytest.raises(SystemExit) as excinfo:
            loader.exec_module(module)
    finally:
        sys.modules.pop(module_name, None)

    message = str(excinfo.value)
    assert "library/cli_utils.py:201" in message
    assert "merge conflict" in message


@pytest.mark.unit
def test_run_stage__inserts_default_document_subcommand(monkeypatch):
    """Document stage should default to the ``all`` subcommand when omitted."""

    from scripts import get_data as cli

    stage = cli.Stage("document", "get_document_data.py")
    tokens = (
        "--log-level",
        "INFO",
        "--limit",
        "5",
        "--input-dir",
        "data\\input",
        "--output-dir",
        "data\\output",
    )
    forward = cli.ForwardArgs(tokens=tokens, extras_start=2, extra_len=6)

    captured: dict[str, list[str]] = {}

    class _Result:
        returncode = 0

    def _fake_run(command, *, check=False, env=None):  # type: ignore[override]
        captured["command"] = command
        return _Result()

    monkeypatch.setattr(cli.subprocess, "run", _fake_run)

    cli.run_stage(stage, forward)

    command = captured["command"]
    assert "all" in command[2:], "expected default 'all' subcommand to be forwarded"
    assert command.index("all") < command.index("--limit")


@pytest.mark.unit
def test_count_output_files__custom_output_directory(tmp_path):
    """Forwarded ``--output-dir`` should be honoured when counting artefacts."""

    from scripts import get_data as cli

    custom_output = tmp_path / "exports"
    custom_output.mkdir()
    (custom_output / "result.csv").write_text("id\n1\n", encoding="utf-8")
    (custom_output / "notes.txt").write_text("skip", encoding="utf-8")

    tokens = ("--log-level", "INFO", f"--output-dir={custom_output}")
    forward = cli.ForwardArgs(tokens=tokens, extras_start=0, extra_len=0)

    resolved = forward.output_dir
    assert resolved == custom_output.resolve()

    count = cli.count_output_files(resolved)
    assert count == 1
