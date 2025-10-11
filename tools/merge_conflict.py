"""Utilities to detect unresolved merge conflict markers."""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Final

__all__ = [
    "MERGE_CONFLICT_PATTERNS",
    "has_merge_conflict_markers",
]


# fmt: off
MERGE_CONFLICT_PATTERNS: Final[tuple[re.Pattern[str], ...]] = (
    re.compile(r"(?m)^[ \t]*<<<<<<<(?=\s)"),
    re.compile(r"(?m)^[ \t]*======= *$"),
    re.compile(r"(?m)^[ \t]*>>>>>>>(?=\s)"),
)
# fmt: on


def has_merge_conflict_markers(text: str | Iterable[str] | bytes) -> bool:
    """Return ``True`` when ``text`` contains merge conflict markers."""

    if isinstance(text, bytes):
        text = text.decode("utf-8", "replace")
    if isinstance(text, Iterable) and not isinstance(text, str):
        text = "".join(str(chunk) for chunk in text)
    if not isinstance(text, str):
        raise TypeError("text must be str, bytes or an iterable producing strings")
    for pattern in MERGE_CONFLICT_PATTERNS:
        if pattern.search(text):
            return True
    return False
