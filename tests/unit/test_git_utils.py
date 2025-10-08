from __future__ import annotations

import importlib
import logging
import subprocess
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


@pytest.mark.unit
def test_git_sha__github_desktop_fallback(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """GitHub Desktop shims should fall back to the bundled Git executable."""

    module = importlib.import_module("library.common.git")

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    repo_root = tmp_path / "repo"
    git_dir = repo_root / ".git"
    repo_root.mkdir()
    git_dir.mkdir()
    (git_dir / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf8")
    ref_dir = git_dir / "refs" / "heads"
    ref_dir.mkdir(parents=True)
    (ref_dir / "main").write_text("deadbeefdeadbeefdeadbeefdeadbeefdeadbeef\n", encoding="utf8")

    desktop_root = tmp_path / "GitHubDesktop"
    stub = desktop_root / "git.exe"
    actual_git = (
        desktop_root
        / "app-3.3.0"
        / "resources"
        / "app"
        / "git"
        / "cmd"
        / "git.exe"
    )
    desktop_root.mkdir()
    actual_git.parent.mkdir(parents=True, exist_ok=True)
    stub.touch()
    actual_git.touch()

    commands: list[list[str]] = []

    def fake_run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(cmd)
        if cmd[0] == str(stub):
            raise subprocess.CalledProcessError(
                returncode=4294967295,
                cmd=cmd,
                stderr="stub failed",
            )
        if cmd[0] == str(actual_git):
            return subprocess.CompletedProcess(cmd, 0, stdout="deadbeef\n", stderr="")
        raise AssertionError(f"Unexpected command: {cmd!r}")

    monkeypatch.setattr(module, "_repo_root", lambda: repo_root)
    monkeypatch.setattr(module.shutil, "which", lambda _: str(stub))
    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.delenv("GIT_SHA", raising=False)

    sha = module._git_sha()

    assert sha == "deadbeef"
    assert [cmd[0] for cmd in commands] == [str(stub), str(actual_git)]

    cache_clear()


@pytest.mark.unit
def test_git_sha__github_desktop_mingw64_fallback(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Modern GitHub Desktop builds bundle Git under ``mingw64/bin``."""

    module = importlib.import_module("library.common.git")

    cache_clear = getattr(module._git_sha, "cache_clear", None)
    if cache_clear is None:  # pragma: no cover - defensive guard
        pytest.skip("_git_sha missing cache_clear helper")
    cache_clear()

    repo_root = tmp_path / "repo"
    git_dir = repo_root / ".git"
    repo_root.mkdir()
    git_dir.mkdir()
    (git_dir / "HEAD").write_text("ref: refs/heads/main\n", encoding="utf8")
    ref_dir = git_dir / "refs" / "heads"
    ref_dir.mkdir(parents=True)
    (ref_dir / "main").write_text("deadbeefdeadbeefdeadbeefdeadbeefdeadbeef\n", encoding="utf8")

    desktop_root = tmp_path / "GitHubDesktop"
    stub = desktop_root / "git.exe"
    cmd_git = (
        desktop_root
        / "app-3.4.0"
        / "resources"
        / "app"
        / "git"
        / "cmd"
        / "git.exe"
    )
    mingw_git = (
        desktop_root
        / "app-3.4.0"
        / "resources"
        / "app"
        / "git"
        / "mingw64"
        / "bin"
        / "git.exe"
    )
    desktop_root.mkdir()
    mingw_git.parent.mkdir(parents=True, exist_ok=True)
    (cmd_git.parent).mkdir(parents=True, exist_ok=True)
    stub.touch()
    cmd_git.touch()
    mingw_git.touch()

    commands: list[list[str]] = []

    def fake_run(cmd: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        commands.append(cmd)
        if cmd[0] == str(stub):
            raise subprocess.CalledProcessError(
                returncode=4294967295,
                cmd=cmd,
                stderr="stub failed",
            )
        if cmd[0] == str(cmd_git):
            raise subprocess.CalledProcessError(
                returncode=1,
                cmd=cmd,
                stderr="cmd shim missing env",
            )
        if cmd[0] == str(mingw_git):
            return subprocess.CompletedProcess(cmd, 0, stdout="cafebabe\n", stderr="")
        raise AssertionError(f"Unexpected command: {cmd!r}")

    monkeypatch.setattr(module, "_repo_root", lambda: repo_root)
    monkeypatch.setattr(module.shutil, "which", lambda _: str(stub))
    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.delenv("GIT_SHA", raising=False)

    sha = module._git_sha()

    assert sha == "cafebabe"
    assert [cmd[0] for cmd in commands] == [str(stub), str(cmd_git), str(mingw_git)]

    cache_clear()
