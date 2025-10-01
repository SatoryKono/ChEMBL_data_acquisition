"""Atomic file writing helpers with cross-platform file locking."""

from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

import portalocker

_LOCK_TIMEOUT_SECONDS = 10.0


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

    with portalocker.Lock(
        str(lock_path),
        mode="w",
        timeout=lock_timeout,
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
