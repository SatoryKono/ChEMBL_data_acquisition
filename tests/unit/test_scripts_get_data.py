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
def test_ensure_base_path_env__injects_default_when_missing():
    """The orchestrator should provide a stable base path to subprocesses."""

    from scripts import get_data as cli

    env: dict[str, str] = {}
    cli._ensure_base_path_env([], env)

    assert env[cli._BASE_PATH_ENV_VAR] == str(cli.PROJECT_ROOT.resolve())


@pytest.mark.unit
def test_ensure_base_path_env__resolves_relative_cli_value(tmp_path, monkeypatch):
    """CLI ``--base-path`` arguments must propagate as absolute paths."""

    from scripts import get_data as cli

    monkeypatch.chdir(tmp_path)
    provided = Path("custom-base")
    env: dict[str, str] = {}

    cli._ensure_base_path_env(["--base-path", str(provided)], env)

    assert env[cli._BASE_PATH_ENV_VAR] == str((tmp_path / provided).resolve())


@pytest.mark.unit
def test_ensure_base_path_env__keeps_preexisting_value():
    """Existing ``CHEMBL_DA_BASE_PATH`` values must remain untouched."""

    from scripts import get_data as cli

    existing = "/tmp/chembl"
    env = {cli._BASE_PATH_ENV_VAR: existing}

    cli._ensure_base_path_env([], env)

    assert env[cli._BASE_PATH_ENV_VAR] == existing


@pytest.mark.unit
def test_cleanup_intermediate_files__removes_known_patterns(tmp_path: Path) -> None:
    """Temporary artefacts should be deleted while canonical files remain."""

    from scripts import get_data as cli

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
def test_cleanup_intermediate_files__skips_missing_directory(tmp_path: Path) -> None:
    """Cleanup should be a no-op when the output directory is absent."""

    from scripts import get_data as cli

    missing_dir = tmp_path / "missing"
    removed = cli.cleanup_intermediate_files(missing_dir)

    assert removed == 0


@pytest.mark.unit
def test_should_run_cleanup__respects_user_flags() -> None:
    """Diagnostic flags should disable the orchestrator cleanup."""

    from scripts import get_data as cli

    base = ("--config", "cfg.yaml", "--output-dir", "data/output")
    forward = cli.ForwardArgs(tokens=base + ("--debug",), extras_start=0, extra_len=len(base) + 1)
    assert cli._should_run_cleanup(forward) is False

    forward = cli.ForwardArgs(tokens=base, extras_start=0, extra_len=len(base))
    assert cli._should_run_cleanup(forward) is True
