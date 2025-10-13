from __future__ import annotations

import importlib
import importlib.util
import logging
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
def test_build_forward_args__respects_equals_style_output_dir(tmp_path):
    """Explicit ``--output-dir=...`` values should not be overridden by defaults."""

    from scripts import get_data as cli

    custom_output = tmp_path / "out"
    args, extras = cli.parse_args([f"--output-dir={custom_output}"])

    forward = cli.build_forward_args(args, extras)

    tokens = forward.tokens
    assert f"--output-dir={custom_output}" in tokens
    assert "--output-dir" not in tokens, "default output-dir flag must not be duplicated"
    assert str(cli.DEFAULT_OUTPUT_DIR) not in tokens


@pytest.mark.unit
def test_main__skip_stage_and_dry_run(monkeypatch, tmp_path):
    """Dry-run skips outputs and does not emit warnings for skipped stages."""

    from scripts import get_data as cli

    base_path = tmp_path / "workspace"
    input_dir = base_path / "data" / "input"
    output_dir = base_path / "data" / "output"
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)

    log_target = tmp_path / "logs" / "orchestrator.log"
    log_target.parent.mkdir(parents=True, exist_ok=True)

    root_logger = logging.getLogger()
    original_level = root_logger.level

    records: list[logging.LogRecord] = []

    class _CaptureHandler(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:  # pragma: no cover - simple sink
            records.append(record)

    capture_handler = _CaptureHandler(level=logging.INFO)
    root_logger.addHandler(capture_handler)

    def _fake_configure_logging(level_name: str | None) -> Path:
        root_logger.setLevel(getattr(logging, level_name or "INFO"))
        return log_target

    monkeypatch.setattr(cli, "configure_logging", _fake_configure_logging)

    invoked: list[str] = []

    def _fake_run_stage(stage, forward_args):  # type: ignore[override]
        invoked.append(stage.name)
        return 0.0

    monkeypatch.setattr(cli, "run_stage", _fake_run_stage)

    config_path = Path(__file__).resolve().parents[2] / "config" / "config.yaml"
    argv = [
        "--skip",
        "activity",
        "--config",
        str(config_path),
        "--dry-run",
        "--base-path",
        str(base_path),
    ]

    try:
        status = cli.main(argv)
    finally:
        root_logger.setLevel(original_level)
        root_logger.removeHandler(capture_handler)

    assert status == 0
    assert "activity" not in invoked
    assert any("Активные этапы (dry-run)" in record.getMessage() for record in records)
    assert any("Dry-run: ожидаемые файлы" in record.getMessage() for record in records)
    assert any(
        "Dry-run: проверка выходных файлов пропущена." in record.getMessage()
        for record in records
    )
    assert not any(record.levelno >= logging.WARNING for record in records)
