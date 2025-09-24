from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any

import pytest

import library.git_utils as git_utils


def _prepare_head(git_dir: Path, commit: str, ref: str = "refs/heads/main") -> None:
    head = git_dir / "HEAD"
    head.parent.mkdir(parents=True, exist_ok=True)
    head.write_text(f"ref: {ref}\n", encoding="utf8")
    ref_path = git_dir / Path(ref)
    ref_path.parent.mkdir(parents=True, exist_ok=True)
    ref_path.write_text(f"{commit}\n", encoding="utf8")


def _configure_failing_git(
    monkeypatch: pytest.MonkeyPatch,
    repo_root: Path,
    *,
    returncode: int = 1,
    stderr: str = "fatal: simulated error",
    stdout: str | None = None,
) -> None:
    def fail(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise subprocess.CalledProcessError(
            returncode=returncode,
            cmd=["git", "rev-parse", "HEAD"],
            stderr=stderr,
            output=stdout,
        )

    monkeypatch.setattr(git_utils, "_repo_root", lambda: repo_root)
    monkeypatch.setattr(git_utils.subprocess, "run", fail)
    monkeypatch.setattr(git_utils.shutil, "which", lambda _: "git")


def _capture_logs(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, dict[str, object]]]:
    records: list[tuple[str, dict[str, object]]] = []

    def info(event: str, *args: object, **kwargs: object) -> None:
        records.append((event, dict(kwargs)))

    monkeypatch.setattr(git_utils.logger, "info", info)
    monkeypatch.setattr(git_utils.logger, "warning", lambda *a, **k: None)
    return records


def _assert_subprocess_payload(
    records: list[tuple[str, dict[str, object]]],
    *,
    expected_returncode: int,
    expected_returncode_raw: int | None = None,
    expected_stderr: str | None = "fatal: simulated error",
) -> None:
    payloads: list[dict[str, Any]] = [
        rec for event, rec in records if event == "git_sha_fallback"
    ]
    assert payloads
    raw = expected_returncode if expected_returncode_raw is None else expected_returncode_raw
    for payload in payloads:
        assert payload.get("error")
        assert payload.get("error_returncode") == expected_returncode
        assert payload.get("error_returncode_raw") == raw
        assert payload.get("error_cmd") == ["git", "rev-parse", "HEAD"]
        if expected_stderr is None:
            assert "error_stderr" not in payload
        else:
            assert payload.get("error_stderr") == expected_stderr


def test_git_sha_fallback_reads_head(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    commit = "1" * 40
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    _prepare_head(git_dir, commit)
    _configure_failing_git(monkeypatch, tmp_path)
    records = _capture_logs(monkeypatch)
    git_utils._git_sha.cache_clear()

    sha = git_utils._git_sha()

    assert sha == commit
    assert any(event == "git_sha_fallback" for event, _ in records)
    assert any(
        rec.get("reason") == "subprocess_error"
        for event, rec in records
        if event == "git_sha_fallback"
    )
    _assert_subprocess_payload(records, expected_returncode=1)


def test_git_sha_fallback_supports_gitdir_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    commit = "2" * 40
    storage = tmp_path / "git_storage"
    git_dir = storage / "worktree"
    git_dir.mkdir(parents=True)
    _prepare_head(git_dir, commit)
    (tmp_path / ".git").write_text("gitdir: git_storage/worktree\n", encoding="utf8")

    _configure_failing_git(monkeypatch, tmp_path)
    records = _capture_logs(monkeypatch)
    git_utils._git_sha.cache_clear()

    sha = git_utils._git_sha()

    assert sha == commit
    assert any(event == "git_sha_fallback" for event, _ in records)
    assert any(
        rec.get("reason") == "subprocess_error"
        for event, rec in records
        if event == "git_sha_fallback"
    )
    _assert_subprocess_payload(records, expected_returncode=1)


def test_git_sha_fallback_reads_packed_refs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    commit = "3" * 40
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    head = git_dir / "HEAD"
    head.write_text("ref: refs/heads/main\n", encoding="utf8")
    packed = git_dir / "packed-refs"
    packed.write_text(f"# comment\n{commit} refs/heads/main\n", encoding="utf8")

    _configure_failing_git(monkeypatch, tmp_path)
    records = _capture_logs(monkeypatch)
    git_utils._git_sha.cache_clear()

    sha = git_utils._git_sha()

    assert sha == commit
    assert any(event == "git_sha_fallback" for event, _ in records)
    assert any(
        rec.get("reason") == "subprocess_error"
        for event, rec in records
        if event == "git_sha_fallback"
    )
    _assert_subprocess_payload(records, expected_returncode=1)


def test_git_sha_fallback_normalises_large_returncode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commit = "5" * 40
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    _prepare_head(git_dir, commit)

    raw_code = 0xFFFFFFFF
    _configure_failing_git(monkeypatch, tmp_path, returncode=raw_code)
    records = _capture_logs(monkeypatch)
    git_utils._git_sha.cache_clear()

    sha = git_utils._git_sha()

    assert sha == commit
    assert any(event == "git_sha_fallback" for event, _ in records)
    _assert_subprocess_payload(
        records,
        expected_returncode=-1,
        expected_returncode_raw=raw_code,
    )


def test_git_sha_fallback_when_git_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    commit = "4" * 40
    git_dir = tmp_path / ".git"
    git_dir.mkdir()
    _prepare_head(git_dir, commit)
    monkeypatch.setattr(git_utils, "_repo_root", lambda: tmp_path)
    monkeypatch.setattr(git_utils.shutil, "which", lambda _: None)
    records = _capture_logs(monkeypatch)
    git_utils._git_sha.cache_clear()

    sha = git_utils._git_sha()

    assert sha == commit
    assert any(event == "git_sha_fallback" for event, _ in records)
    assert any(
        rec.get("reason") == "missing_executable"
        for event, rec in records
        if event == "git_sha_fallback"
    )
