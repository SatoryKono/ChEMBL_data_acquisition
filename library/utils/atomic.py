"""Atomic file writing helpers with cross-platform file locking."""

from __future__ import annotations

import errno
import os
import tempfile
import time
from collections.abc import Iterator
from contextlib import AbstractContextManager, contextmanager
from pathlib import Path
from types import ModuleType, TracebackType
from typing import IO, TYPE_CHECKING, Any, Protocol, TextIO, cast

__all__ = ["open_atomic", "robust_replace"]

_REPLACE_RETRYABLE_ERRNOS = {errno.EACCES, errno.EPERM}
_REPLACE_RETRYABLE_WINERRORS = {5, 32, 33}


class _PortalockerLockFallback(AbstractContextManager[Any]):  # pragma: no cover - helper
    ...


if TYPE_CHECKING:  # pragma: no cover - typing helpers only
    try:
        from portalocker import (
            Lock as _PortalockerLock,  # type: ignore[import-not-found, unused-ignore]  # pragma: no cover
        )
    except ModuleNotFoundError:  # pragma: no cover - typing fallback
        _PortalockerLock = _PortalockerLockFallback
else:  # pragma: no cover - runtime helper
    _PortalockerLock = _PortalockerLockFallback


class _PortalockerModule(Protocol):
    def Lock(
        self,
        path: str,
        *,
        mode: str = ...,
        timeout: float = ...,
    ) -> _PortalockerLock: ...

try:  # pragma: no cover - optional dependency
    import portalocker as _portalocker  # type: ignore[import-not-found, unused-ignore]
except ModuleNotFoundError:  # pragma: no cover - fallback path
    portalocker: _PortalockerModule | None = None
else:
    portalocker = cast(_PortalockerModule, _portalocker)

if TYPE_CHECKING:  # pragma: no cover - typing helpers only

    class _MSVCRTModule(Protocol):
        LK_NBLCK: int
        LK_UNLCK: int

        def locking(self, fd: int, mode: int, size: int) -> None: ...

    class _FcntlModule(Protocol):
        LOCK_EX: int
        LOCK_NB: int
        LOCK_UN: int

        def flock(self, fd: int, op: int) -> None: ...

else:  # pragma: no cover - runtime branch
    _MSVCRTModule = ModuleType
    _FcntlModule = ModuleType

if os.name == "nt":  # pragma: no cover - platform dependent
    import msvcrt as _msvcrt

    msvcrt: _MSVCRTModule | None = cast("_MSVCRTModule", _msvcrt)
    fcntl: _FcntlModule | None = None
else:  # pragma: no cover - platform dependent
    msvcrt = None
    import fcntl as _fcntl

    fcntl = cast("_FcntlModule", _fcntl)

_LOCK_TIMEOUT_SECONDS = 10.0
_LOCK_POLL_INTERVAL_SECONDS = 0.1


class FileLock:
    """Lightweight file lock implementation using OS primitives."""

    def __init__(self, path: Path, *, timeout: float, poll_interval: float) -> None:
        self._path = path
        self._timeout = timeout
        self._poll_interval = poll_interval
        self._handle: TextIO | None = None

    def __enter__(self) -> FileLock:
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

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None:
        if self._handle is None:
            return

        _unlock(self._handle)
        self._handle.close()
        self._handle = None


def _try_lock(handle: TextIO) -> bool:
    """Attempt to acquire an exclusive lock on *handle* without blocking."""

    try:
        if os.name == "nt":  # pragma: no cover - platform-specific behaviour
            if msvcrt is None:  # pragma: no cover - defensive branch
                raise RuntimeError("msvcrt is unavailable on this platform")
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:  # pragma: no cover - platform-specific behaviour
            if fcntl is None:  # pragma: no cover - defensive branch
                raise RuntimeError("fcntl is unavailable on this platform")
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
        if msvcrt is None:  # pragma: no cover - defensive branch
            raise RuntimeError("msvcrt is unavailable on this platform")
        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    else:  # pragma: no cover - platform-specific behaviour
        if fcntl is None:  # pragma: no cover - defensive branch
            raise RuntimeError("fcntl is unavailable on this platform")
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
) -> Iterator[IO[Any]]:
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
            tmp_file: IO[Any]
            tmp_file = os.fdopen(fd, mode, encoding=encoding, newline=newline)
            with tmp_file:
                yield tmp_file
        except Exception:
            raise
        else:
            robust_replace(tmp_name, path)
        finally:
            try:
                os.unlink(tmp_name)
            except FileNotFoundError:
                pass


def robust_replace(
    src: str | os.PathLike[str],
    dst: str | os.PathLike[str],
    *,
    attempts: int = 5,
    delay: float = 0.1,
    backoff: float = 2.0,
) -> None:
    """Replace *dst* with *src* retrying transient permission errors."""

    if attempts <= 0:
        raise ValueError("attempts must be a positive integer")
    if delay < 0:
        raise ValueError("delay must be non-negative")
    if backoff < 1:
        raise ValueError("backoff must be greater than or equal to 1")

    attempt = 0
    sleep_delay = delay
    while True:
        attempt += 1
        try:
            os.replace(src, dst)
            return
        except OSError as exc:
            if attempt >= attempts or not _should_retry_replace(exc):
                raise

            time.sleep(sleep_delay)
            sleep_delay *= backoff


def _should_retry_replace(exc: OSError) -> bool:
    """Return ``True`` if ``exc`` represents a transient replace failure."""

    err_no = exc.errno
    win_err = getattr(exc, "winerror", None)
    return (err_no in _REPLACE_RETRYABLE_ERRNOS) or (
        win_err in _REPLACE_RETRYABLE_WINERRORS
    )


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
            time.sleep(_LOCK_POLL_INTERVAL_SECONDS)
        except OSError as exc:  # pragma: no cover - platform dependent
            if os.name == "nt":
                win_err = getattr(exc, "winerror", None)
                if win_err not in {32, 33}:  # sharing violation / lock violation
                    raise
            else:
                if exc.errno not in {errno.EACCES, errno.EAGAIN}:
                    raise
            if timeout is not None and time.monotonic() - start_time > timeout:
                raise TimeoutError("timed out acquiring file lock") from exc
            time.sleep(_LOCK_POLL_INTERVAL_SECONDS)
        else:
            break


def _lock_file(lock_file: TextIO) -> None:
    """Lock the file in a platform-specific manner."""

    if os.name == "nt":  # pragma: no branch - platform specific
        if msvcrt is None:  # pragma: no cover - defensive branch
            raise RuntimeError("msvcrt is unavailable on this platform")
        msvcrt.locking(lock_file.fileno(), msvcrt.LK_NBLCK, 1)
    else:
        if fcntl is None:  # pragma: no cover - defensive branch
            raise RuntimeError("fcntl is unavailable on this platform")
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def _release_lock(lock_file: TextIO) -> None:
    """Release the platform-specific file lock."""

    if os.name == "nt":  # pragma: no branch - platform specific
        if msvcrt is None:  # pragma: no cover - defensive branch
            raise RuntimeError("msvcrt is unavailable on this platform")
        msvcrt.locking(lock_file.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        if fcntl is None:  # pragma: no cover - defensive branch
            raise RuntimeError("fcntl is unavailable on this platform")
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


def _prepare_lock_file(lock_file: TextIO) -> None:
    """Ensure the lock file has a byte to lock and reset the pointer."""

    lock_file.seek(0, os.SEEK_END)
    if lock_file.tell() == 0:
        lock_file.write("0")
        lock_file.flush()
    lock_file.seek(0)
