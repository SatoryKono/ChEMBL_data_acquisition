"""Unit tests for :mod:`library.utils.atomic`."""

from __future__ import annotations

import errno
from pathlib import Path

import pytest

from library.utils import atomic


@pytest.mark.unit
def test_robust_replace__retries_on_transient_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    src = tmp_path / "source.txt"
    dst = tmp_path / "destination.txt"
    src.write_text("new", encoding="utf-8")
    dst.write_text("old", encoding="utf-8")

    real_replace = atomic.os.replace
    calls: list[tuple[str, str]] = []

    def flaky_replace(src_path: str | Path, dst_path: str | Path) -> None:
        calls.append((str(src_path), str(dst_path)))
        if len(calls) < 3:
            raise PermissionError(errno.EACCES, "Access is denied", str(dst_path))
        real_replace(str(src_path), str(dst_path))

    monkeypatch.setattr(atomic.os, "replace", flaky_replace)
    monkeypatch.setattr(atomic.time, "sleep", lambda _seconds: None)

    atomic.robust_replace(src, dst, attempts=5, delay=0.01, backoff=1.0)

    assert len(calls) == 3
    assert dst.read_text(encoding="utf-8") == "new"
    assert not src.exists()


@pytest.mark.unit
def test_robust_replace__raises_after_exhausting_attempts(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    src = tmp_path / "stuck.txt"
    dst = tmp_path / "blocked.txt"
    src.write_text("payload", encoding="utf-8")

    calls: list[tuple[str, str]] = []

    def always_fail(src_path: str | Path, dst_path: str | Path) -> None:
        calls.append((str(src_path), str(dst_path)))
        raise PermissionError(errno.EACCES, "Access is denied", str(dst_path))

    monkeypatch.setattr(atomic.os, "replace", always_fail)
    monkeypatch.setattr(atomic.time, "sleep", lambda _seconds: None)

    with pytest.raises(PermissionError):
        atomic.robust_replace(src, dst, attempts=2, delay=0.01, backoff=1.0)

    assert len(calls) == 2
    assert src.exists()
    assert not dst.exists()
