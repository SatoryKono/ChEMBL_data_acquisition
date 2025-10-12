from __future__ import annotations

import importlib
import importlib.util
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path
from types import SimpleNamespace

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
    sys.modules[module_name] = module

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
def test_build_forward_args__injects_project_paths(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from scripts import get_data

    base_path = tmp_path / "data-root"
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))

    args = SimpleNamespace(limit=None, log_level=None, config=None)
    forwarded = get_data.build_forward_args(args, [])

    tokens = list(forwarded.tokens)
    assert "--base-path" in tokens
    base_index = tokens.index("--base-path")
    assert Path(tokens[base_index + 1]) == base_path

    expected_tail = [
        "--input-dir",
        get_data.DEFAULT_INPUT_DIR,
        "--output-dir",
        get_data.DEFAULT_OUTPUT_DIR,
    ]
    assert tokens[-len(expected_tail) :] == expected_tail


@pytest.mark.unit
def test_build_forward_args__keeps_user_supplied_paths() -> None:
    from scripts import get_data

    args = SimpleNamespace(limit=None, log_level=None, config=None)
    extra = [
        "--base-path",
        "/custom/base",
        "--input-dir",
        "incoming",
        "--output-dir",
        "processed",
    ]

    forwarded = get_data.build_forward_args(args, extra)

    assert forwarded.tokens.count("--base-path") == 1
    assert forwarded.tokens.count("--input-dir") == 1
    assert forwarded.tokens.count("--output-dir") == 1
    assert list(forwarded.tokens[: len(extra)]) == extra


@pytest.mark.unit
def test_parse_args__coerces_positional_limit() -> None:
    from scripts import get_data

    args, unknown = get_data.parse_args(["10"])

    assert args.limit == 10
    assert unknown == []


@pytest.mark.unit
def test_build_forward_args__target_command_precedes_limit() -> None:
    from scripts import get_data

    args = SimpleNamespace(limit=7, log_level=None, config=None)
    forwarded = get_data.build_forward_args(args, [])

    stage_args = forwarded.with_default_subcommand(
        "all", choices=get_data.TARGET_SUBCOMMANDS
    )

    assert stage_args.count("all") == 1
    assert "--limit" in stage_args
    assert stage_args.index("--limit") > stage_args.index("all")


@pytest.mark.unit
def test_parse_args__positional_limit_conflict(capsys: pytest.CaptureFixture[str]) -> None:
    from scripts import get_data

    with pytest.raises(SystemExit) as excinfo:
        get_data.parse_args(["--limit", "5", "10"])

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Ограничение уже указано" in captured.err


@pytest.mark.unit
def test_parse_args__positional_limit_with_extra_tokens(capsys: pytest.CaptureFixture[str]) -> None:
    from scripts import get_data

    with pytest.raises(SystemExit) as excinfo:
        get_data.parse_args(["10", "chembl"])

    assert excinfo.value.code == 2
    captured = capsys.readouterr()
    assert "Неизвестный позиционный аргумент" in captured.err


@pytest.mark.unit
def test_run_stage__forces_pubchem_env_override(monkeypatch, tmp_path) -> None:
    from scripts import get_data

    stage = get_data.Stage("testitem", "get_testitem_data.py")
    script_path = tmp_path / stage.script
    script_path.write_text("print('ok')", encoding="utf-8")

    monkeypatch.setattr(get_data, "SCRIPTS_DIR", tmp_path, raising=False)
    monkeypatch.delenv("CHEMBL_DA_PUBCHEM_ENABLE", raising=False)

    captured_env: dict[str, str] = {}

    class _Result:
        returncode = 0

    def _fake_run(command, check, env):  # type: ignore[override]
        captured_env.update(env)
        return _Result()

    config = SimpleNamespace(
        sources=SimpleNamespace(pubchem=SimpleNamespace(enable=False))
    )

    def _fake_load_config(path, base_path=None):  # type: ignore[override]
        return config

    monkeypatch.setattr(get_data.subprocess, "run", _fake_run)
    monkeypatch.setattr(get_data, "load_config", _fake_load_config, raising=False)

    get_data.run_stage(stage, ["--config", str(tmp_path / "config.yaml")])

    assert captured_env.get("CHEMBL_DA_PUBCHEM_ENABLE") == "true"


@pytest.mark.unit
def test_run_stage__respects_enabled_config(monkeypatch, tmp_path) -> None:
    from scripts import get_data

    stage = get_data.Stage("testitem", "get_testitem_data.py")
    script_path = tmp_path / stage.script
    script_path.write_text("print('ok')", encoding="utf-8")

    monkeypatch.setattr(get_data, "SCRIPTS_DIR", tmp_path, raising=False)
    monkeypatch.delenv("CHEMBL_DA_PUBCHEM_ENABLE", raising=False)

    class _Result:
        returncode = 0

    captured_env: dict[str, str] = {}

    def _fake_run(command, check, env):  # type: ignore[override]
        captured_env.update(env)
        return _Result()

    config = SimpleNamespace(
        sources=SimpleNamespace(pubchem=SimpleNamespace(enable=True))
    )

    def _fake_load_config(path, base_path=None):  # type: ignore[override]
        return config

    monkeypatch.setattr(get_data.subprocess, "run", _fake_run)
    monkeypatch.setattr(get_data, "load_config", _fake_load_config, raising=False)

    get_data.run_stage(stage, ["--config", str(tmp_path / "config.yaml")])

    assert "CHEMBL_DA_PUBCHEM_ENABLE" not in captured_env


@pytest.mark.unit
def test_run_stage__populates_base_path_env(monkeypatch, tmp_path) -> None:
    from scripts import get_data

    stage = get_data.Stage("document", "get_document_data.py")
    script_path = tmp_path / stage.script
    script_path.write_text("print('ok')", encoding="utf-8")

    monkeypatch.setattr(get_data, "SCRIPTS_DIR", tmp_path, raising=False)
    monkeypatch.chdir(tmp_path)

    captured_env: dict[str, str] = {}

    class _Result:
        returncode = 0

    def _fake_run(command, check, env):  # type: ignore[override]
        captured_env.update(env)
        return _Result()

    monkeypatch.setattr(get_data.subprocess, "run", _fake_run)

    get_data.run_stage(stage, ["--base-path", "workspace"])

    expected_base = (tmp_path / "workspace").resolve()
    assert captured_env.get("CHEMBL_DA_BASE_PATH") == str(expected_base)


@pytest.mark.unit
def test_run_stage__preserves_existing_base_path_env(monkeypatch, tmp_path) -> None:
    from scripts import get_data

    stage = get_data.Stage("document", "get_document_data.py")
    script_path = tmp_path / stage.script
    script_path.write_text("print('ok')", encoding="utf-8")

    monkeypatch.setattr(get_data, "SCRIPTS_DIR", tmp_path, raising=False)
    monkeypatch.chdir(tmp_path)

    captured_env: dict[str, str] = {}

    class _Result:
        returncode = 0

    def _fake_run(command, check, env):  # type: ignore[override]
        captured_env.update(env)
        return _Result()

    monkeypatch.setattr(get_data.subprocess, "run", _fake_run)
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(tmp_path / "preset"))

    get_data.run_stage(stage, ["--base-path", "workspace"])

    assert captured_env.get("CHEMBL_DA_BASE_PATH") == str(tmp_path / "preset")
