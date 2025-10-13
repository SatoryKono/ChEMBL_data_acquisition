from __future__ import annotations

import importlib
import importlib.util
import logging
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

import pytest


@pytest.fixture()
def cli(monkeypatch: pytest.MonkeyPatch):
    """Import ``scripts.get_data`` with the heavy CLI module stubbed out."""

    import types

    module_name = "library.cli.commands.get_data"
    if module_name not in sys.modules:
        stub = types.ModuleType(module_name)
        monkeypatch.setitem(sys.modules, module_name, stub)
    module = importlib.import_module("scripts.get_data")
    return module


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
def test_run_stage__inserts_default_document_subcommand(
    monkeypatch: pytest.MonkeyPatch, cli
):
    """Document stage should default to the ``all`` subcommand when omitted."""

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
def test_build_forward_args__respects_equals_style_output_dir(
    tmp_path: Path, cli
):
    """Explicit ``--output-dir=...`` values should not be overridden by defaults."""

    custom_output = tmp_path / "out"
    args, extras = cli.parse_args([f"--output-dir={custom_output}"])

    forward = cli.build_forward_args(args, extras)

    tokens = forward.tokens
    assert f"--output-dir={custom_output}" in tokens
    assert "--output-dir" not in tokens, "default output-dir flag must not be duplicated"
    assert str(cli.DEFAULT_OUTPUT_DIR) not in tokens


@pytest.mark.unit
def test_ensure_base_path_env__injects_default_when_missing(cli) -> None:
    """The orchestrator should provide a stable base path to subprocesses."""

    env: dict[str, str] = {}
    cli._ensure_base_path_env([], env)

    assert env[cli._BASE_PATH_ENV_VAR] == str(cli.PROJECT_ROOT.resolve())


@pytest.mark.unit
def test_ensure_base_path_env__resolves_relative_cli_value(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cli
):
    """CLI ``--base-path`` arguments must propagate as absolute paths."""

    monkeypatch.chdir(tmp_path)
    provided = Path("custom-base")
    env: dict[str, str] = {}

    cli._ensure_base_path_env(["--base-path", str(provided)], env)

    assert env[cli._BASE_PATH_ENV_VAR] == str((tmp_path / provided).resolve())


@pytest.mark.unit
def test_ensure_base_path_env__keeps_preexisting_value(cli) -> None:
    """Existing ``CHEMBL_DA_BASE_PATH`` values must remain untouched."""

    existing = "/tmp/chembl"
    env = {cli._BASE_PATH_ENV_VAR: existing}

    cli._ensure_base_path_env([], env)

    assert env[cli._BASE_PATH_ENV_VAR] == existing


@pytest.mark.unit
def test_cleanup_intermediate_files__removes_known_patterns(
    tmp_path: Path, cli
) -> None:
    """Temporary artefacts should be deleted while canonical files remain."""

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    removable_files = [
        output_dir / ".output.targets_20240101.csv.tmp",
        output_dir / "output.targets_20240101_intermediate.csv",
        output_dir / "activity_debug.log",
        output_dir / "cache_chunk.pkl",
        output_dir / "pubchem_request.json.lock",
    ]
    for path in removable_files:
        path.write_text("dummy", encoding="utf-8")

    raw_dir = output_dir / "raw"
    raw_dir.mkdir()
    (raw_dir / "temporary.csv").write_text("id\n1\n", encoding="utf-8")

    canonical = output_dir / "output.targets_20240101.csv"
    canonical.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    removed = cli.cleanup_intermediate_files(output_dir)

    assert removed == len(removable_files) + 1  # +1 for the raw directory
    assert canonical.exists()
    assert not raw_dir.exists()
    assert all(not path.exists() for path in removable_files)


@pytest.mark.unit
def test_cleanup_intermediate_files__skips_missing_directory(
    tmp_path: Path, cli
) -> None:
    """Cleanup should be a no-op when the output directory is absent."""

    missing_dir = tmp_path / "missing"
    removed = cli.cleanup_intermediate_files(missing_dir)

    assert removed == 0


@pytest.mark.unit
def test_should_run_cleanup__respects_user_flags(cli) -> None:
    """Diagnostic flags should disable the orchestrator cleanup."""

    base = ("--config", "cfg.yaml", "--output-dir", "data/output")
    forward = cli.ForwardArgs(tokens=base + ("--debug",), extras_start=0, extra_len=len(base) + 1)
    assert cli._should_run_cleanup(forward) is False

    forward = cli.ForwardArgs(tokens=base, extras_start=0, extra_len=len(base))
    assert cli._should_run_cleanup(forward) is True


@pytest.mark.unit
def test_main_skip_stage__does_not_emit_output_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
    cli,
) -> None:
    """Skipping a stage should not trigger output count warnings."""

    output_dir = tmp_path / "output"
    output_dir.mkdir()
    logs_dir = tmp_path / "logs"

    monkeypatch.setattr(cli, "OUTPUT_DIR", output_dir)
    monkeypatch.setattr(cli, "LOGS_DIR", logs_dir)

    executed: list[str] = []

    def _fake_run_stage(stage: cli.Stage, forward_args: cli.ForwardArgs) -> float:
        executed.append(stage.name)
        artefact = output_dir / f"{stage.name}_{len(executed)}.csv"
        artefact.write_text("id\n1\n", encoding="utf-8")
        return 1.0

    monkeypatch.setattr(cli, "run_stage", _fake_run_stage)

    cleanup_calls: list[Path] = []

    def _fake_cleanup(path: Path) -> int:
        cleanup_calls.append(path)
        return 0

    monkeypatch.setattr(cli, "cleanup_intermediate_files", _fake_cleanup)

    def _fake_configure_logging(level_name: str | None) -> Path:
        logs_dir.mkdir(parents=True, exist_ok=True)
        log_path = logs_dir / "get_data.log"
        log_path.write_text("", encoding="utf-8")
        level = (
            logging.INFO
            if level_name is None
            else logging._nameToLevel.get(level_name, logging.INFO)
        )
        logging.getLogger().setLevel(level)
        return log_path

    monkeypatch.setattr(cli, "configure_logging", _fake_configure_logging)

    caplog.set_level(logging.INFO)

    exit_code = cli.main(["--skip", "activity"])

    assert exit_code == 0
    assert executed == ["testitem", "target", "document", "assay"]
    assert cleanup_calls == [output_dir]
    csv_files = sorted(path.name for path in output_dir.glob("*.csv"))
    assert len(csv_files) == len(executed)

    warning_messages = [
        record.message for record in caplog.records if record.levelno >= logging.WARNING
    ]
    assert not warning_messages
    assert any("Найдено" in record.message for record in caplog.records)
