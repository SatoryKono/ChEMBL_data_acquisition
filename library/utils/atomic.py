"""Atomic file writing helpers with cross-platform file locking."""

from __future__ import annotations

import errno
import os
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

try:  # pragma: no cover - optional dependency
    import portalocker
except ModuleNotFoundError:  # pragma: no cover - fallback path
    portalocker = None

if os.name == "nt":  # pragma: no cover - platform dependent
    import msvcrt
else:  # pragma: no cover - platform dependent
    import errno

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

    with _acquire_lock(lock_path, lock_timeout):

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


@contextmanager
def _acquire_lock(lock_path: Path, timeout: float) -> Iterator[None]:
    """Acquire a filesystem lock with an optional portalocker fallback."""

    if portalocker is not None:
        with portalocker.Lock(
            str(lock_path),
            mode="w",
            timeout=timeout,
        ):
            yield
        return

    start = time.monotonic()
    with open(lock_path, "a+", encoding="utf-8") as lock_file:
        _prepare_lock_file(lock_file)
        _wait_for_lock(lock_file, timeout, start)
        try:
            yield
        finally:
            _release_lock(lock_file)


def _wait_for_lock(lock_file: TextIO, timeout: float, start_time: float) -> None:
    """Attempt to obtain the lock, respecting the timeout."""

    while True:
        try:
            _lock_file(lock_file)
        except BlockingIOError as exc:  # pragma: no cover - dependent on timing
            if timeout is not None and time.monotonic() - start_time > timeout:
                raise TimeoutError("timed out acquiring file lock") from exc
            time.sleep(0.1)
        except OSError as exc:  # pragma: no cover - platform dependent
            if os.name == "nt":
                if exc.winerror not in {32, 33}:  # sharing violation / lock violation
                    raise
            else:
                if exc.errno not in {errno.EACCES, errno.EAGAIN}:
                    raise
            if timeout is not None and time.monotonic() - start_time > timeout:
                raise TimeoutError("timed out acquiring file lock") from exc
            time.sleep(0.1)
        else:
            break


def _lock_file(lock_file: TextIO) -> None:
    """Lock the file in a platform-specific manner."""

    if os.name == "nt":  # pragma: no branch - platform specific
        msvcrt.locking(lock_file.fileno(), msvcrt.LK_NBLCK, 1)
    else:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def _release_lock(lock_file: TextIO) -> None:
    """Release the platform-specific file lock."""

    if os.name == "nt":  # pragma: no branch - platform specific
        msvcrt.locking(lock_file.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _prepare_lock_file(lock_file: TextIO) -> None:
    """Ensure the lock file has a byte to lock and reset the pointer."""

    lock_file.seek(0, os.SEEK_END)
    if lock_file.tell() == 0:
        lock_file.write("0")
        lock_file.flush()
    lock_file.seek(0)
