"""Atomic file writing helpers with cross-platform file locking."""

from __future__ import annotations

import errno
import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

if os.name == "nt":  # pragma: no cover - platform-specific import
    import msvcrt
else:  # pragma: no cover - platform-specific import
    import fcntl

_LOCK_TIMEOUT_SECONDS = 10.0
_LOCK_POLL_INTERVAL_SECONDS = 0.1


class FileLock:
    """Lightweight file lock implementation using OS primitives."""

    def __init__(self, path: Path, *, timeout: float, poll_interval: float) -> None:
        self._path = path
        self._timeout = timeout
        self._poll_interval = poll_interval
        self._handle: TextIO | None = None

    def __enter__(self) -> "FileLock":
        start_time = time.monotonic()
        self._path.parent.mkdir(parents=True, exist_ok=True)
        handle = self._path.open("w+")

        while True:
            if _try_lock(handle):
                self._handle = handle
                return self

            elapsed = time.monotonic() - start_time
            if elapsed >= self._timeout:
                handle.close()
                msg = (
                    f"Could not acquire lock for '{self._path}' within "
                    f"{self._timeout:.2f} seconds"
                )
                raise TimeoutError(msg)

            time.sleep(self._poll_interval)

    def __exit__(self, exc_type, exc, tb) -> None:  # type: ignore[override]
        if self._handle is None:
            return

        _unlock(self._handle)
        self._handle.close()
        self._handle = None


def _try_lock(handle: TextIO) -> bool:
    """Attempt to acquire an exclusive lock on *handle* without blocking."""

    try:
        if os.name == "nt":  # pragma: no cover - platform-specific behaviour
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:  # pragma: no cover - platform-specific behaviour
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError as exc:
        if _is_lock_contended(exc):
            return False
        raise
    except BlockingIOError:
        return False
    return True


def _unlock(handle: TextIO) -> None:
    """Release a previously acquired lock on *handle*."""

    if os.name == "nt":  # pragma: no cover - platform-specific behaviour
        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    else:  # pragma: no cover - platform-specific behaviour
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _is_lock_contended(exc: OSError) -> bool:
    """Return ``True`` when *exc* indicates the lock is already held."""

    err_no = exc.errno
    win_err = getattr(exc, "winerror", None)
    return err_no in {errno.EACCES, errno.EAGAIN} or win_err == 33


@contextmanager
def open_atomic(
    path: Path,
    *,
    mode: str = "w",
    encoding: str | None = None,
    newline: str | None = None,
    lock_timeout: float = _LOCK_TIMEOUT_SECONDS,
) -> Iterator[TextIO]:
    """Open *path* for atomic writing.

    A temporary file is created in the same directory and moved into place with
    :func:`os.replace` after the context manager exits successfully. A
    ``.lock`` sidecar file is used to coordinate concurrent writers.
    """

    if "w" not in mode and "+" not in mode and "a" not in mode:
        msg = "atomic writes require a write-capable mode"
        raise ValueError(msg)

    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    lock_path = path.with_name(f"{path.name}.lock")

    with FileLock(
        lock_path,
        timeout=lock_timeout,
        poll_interval=_LOCK_POLL_INTERVAL_SECONDS,
    ):
        fd, tmp_name = tempfile.mkstemp(
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            text="b" not in mode,
        )
        try:
            with os.fdopen(fd, mode, encoding=encoding, newline=newline) as tmp_file:
                yield tmp_file
        except Exception:
            raise
        else:
            os.replace(tmp_name, path)
        finally:
            try:
                os.unlink(tmp_name)
            except FileNotFoundError:
                pass
