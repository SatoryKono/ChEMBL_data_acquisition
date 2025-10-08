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
def test_git_sha__prefers_head_lookup(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    module_name: str,
) -> None:
    """Reading ``HEAD`` directly should avoid spawning ``git`` subprocesses."""

    module = importlib.import_module(module_name)

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    target_module = importlib.import_module(module._git_sha.__module__)

    repo_root = tmp_path / "project"
    git_dir = repo_root / ".git"
    head_ref = git_dir / "refs" / "heads"
    head_ref.mkdir(parents=True)

    expected_sha = "3f87ce9ea06b4d28efddabd2203fbb1fee117a06"
    (git_dir / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf8")
    (head_ref / "main").write_text(f"{expected_sha}\n", encoding="utf8")

    monkeypatch.setattr(target_module, "_repo_root", lambda: repo_root)

    # ``subprocess.run`` must not be invoked when ``HEAD`` already yields a SHA.
    def _unexpected_run(*_: object, **__: object) -> None:  # pragma: no cover - defensive
        raise AssertionError("subprocess.run should not be called when HEAD is readable")

    monkeypatch.setattr(target_module.subprocess, "run", _unexpected_run)

    sha = module._git_sha()

    assert sha == expected_sha


@pytest.mark.unit
def test_git_sha__shared_singleton() -> None:
    """Both modules should expose the same cached helper instance."""

    common = importlib.import_module("library.common.git")
    legacy = importlib.import_module("library.git_utils")

    assert legacy._git_sha is common._git_sha
