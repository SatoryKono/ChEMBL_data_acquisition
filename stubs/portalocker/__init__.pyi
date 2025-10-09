from __future__ import annotations

import os
from contextlib import AbstractContextManager
from typing import Any

__all__ = ["Lock"]


def Lock(
    file: str | os.PathLike[str] | int,
    *,
    mode: str = ...,
    timeout: float = ...,
    flags: int = ...,
    check_interval: float = ...,
    fail_when_locked: bool = ...,
) -> AbstractContextManager[Any]: ...
