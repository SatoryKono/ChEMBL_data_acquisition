from __future__ import annotations

import importlib
import logging
from pathlib import Path

import pytest


def _expected_repo_root() -> Path:
    current = Path(__file__).resolve().parent
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists() or (candidate / "pyproject.toml").is_file():
            return candidate
    raise AssertionError("Unable to determine repository root for tests")


@pytest.mark.unit
@pytest.mark.parametrize("module_name", ["library.git_utils", "library.common.git"])
def test_repo_root__detects_project_root(module_name: str) -> None:
    """``_repo_root`` should resolve to the actual repository root."""

    module = importlib.import_module(module_name)

    assert module._repo_root() == _expected_repo_root()


@pytest.mark.unit
@pytest.mark.parametrize("module_name", ["library.git_utils", "library.common.git"])
def test_git_sha__missing_git_directory_logged_as_info(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    module_name: str,
) -> None:
    """Missing ``.git`` directories should not emit warning-level logs."""

    module = importlib.import_module(module_name)

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    target_module = importlib.import_module(module._git_sha.__module__)

    monkeypatch.setattr(target_module, "_repo_root", lambda: tmp_path)

    caplog.set_level(logging.DEBUG, logger="chembl")

    sha = module._git_sha()

    assert sha == "UNKNOWN"

    directory_logs = [
        (record.levelno, record.message)
        for record in caplog.records
        if "git_directory_missing" in record.message
    ]

    assert directory_logs, "expected git_directory_missing log entry"
    assert all(level < logging.WARNING for level, _ in directory_logs)


@pytest.mark.unit
@pytest.mark.parametrize("module_name", ["library.git_utils", "library.common.git"])
def test_git_sha__prefers_head_file_over_subprocess(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    module_name: str,
) -> None:
    """Reading ``HEAD`` directly should avoid invoking the Git executable."""

    module = importlib.import_module(module_name)

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    target_module = importlib.import_module(module._git_sha.__module__)

    monkeypatch.setattr(target_module, "_repo_root", lambda: tmp_path)

    git_dir = tmp_path / ".git"
    ref_dir = git_dir / "refs" / "heads"
    ref_dir.mkdir(parents=True)

    sha = "a" * 40
    (git_dir / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf8")
    (ref_dir / "main").write_text(f"{sha}\n", encoding="utf8")

    def _fail_git(*args: object, **kwargs: object) -> None:
        raise AssertionError("git executable should not be invoked")

    monkeypatch.setattr(target_module.shutil, "which", lambda _: "git")
    monkeypatch.setattr(target_module.subprocess, "run", _fail_git)

    caplog.set_level(logging.DEBUG, logger="chembl")

    result = module._git_sha()

    assert result == sha

    head_logs = [
        record for record in caplog.records if "git_sha_head" in record.message
    ]
    assert head_logs, "expected git_sha_head log entry"


@pytest.mark.unit
def test_git_sha__shared_singleton() -> None:
    """Both modules should expose the same cached helper instance."""

    common = importlib.import_module("library.common.git")
    legacy = importlib.import_module("library.git_utils")

    assert legacy._git_sha is common._git_sha
