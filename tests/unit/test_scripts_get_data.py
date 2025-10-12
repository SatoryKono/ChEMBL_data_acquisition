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
def test_build_forward_args__injects_project_paths() -> None:
    from scripts import get_data

    args = SimpleNamespace(limit=None, log_level=None, config=None)
    forwarded = get_data.build_forward_args(args, [])

    expected_tail = [
        "--base-path",
        str(get_data.DATA_DIR),
        "--input-dir",
        get_data.DEFAULT_INPUT_DIR,
        "--output-dir",
        get_data.DEFAULT_OUTPUT_DIR,
    ]
    assert forwarded[-len(expected_tail) :] == expected_tail


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

    assert forwarded.count("--base-path") == 1
    assert forwarded.count("--input-dir") == 1
    assert forwarded.count("--output-dir") == 1
    assert forwarded[: len(extra)] == extra
