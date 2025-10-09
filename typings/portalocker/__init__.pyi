from __future__ import annotations

from contextlib import AbstractContextManager
from types import TracebackType
from typing import Any


class Lock(AbstractContextManager[Any]):
    """Typed stub for :func:`portalocker.Lock`."""

    path: str
    mode: str
    timeout: float

    def __init__(
        self,
        filename: str,
        mode: str = ...,
        timeout: float = ...,
        *,
        flags: int | None = ...,
        check_interval: float | None = ...,
        fail_when_locked: bool | None = ...,
    ) -> None: ...

    def acquire(
        self,
        timeout: float | None = ...,
        check_interval: float | None = ...,
        fail_when_locked: bool | None = ...,
    ) -> bool: ...

    def release(self) -> None: ...

    def __enter__(self) -> Lock: ...

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        tb: TracebackType | None,
    ) -> None: ...


__all__ = ["Lock"]
