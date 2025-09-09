"""Shared utility functions."""

from __future__ import annotations

from typing import Sequence


def pipe_merge(values: Sequence[str | None]) -> str:
    """Return a ``"|"``-joined string of unique, non-empty tokens.

    Parameters
    ----------
    values:
        Sequence of pipe-delimited strings to merge.

    Returns
    -------
    str
        Sorted, unique tokens separated by ``"|"``. Empty inputs yield
        an empty string.
    """

    tokens: set[str] = set()
    for value in values:
        if isinstance(value, str) and value:
            parts = [p.strip() for p in value.split("|") if p.strip()]
            tokens.update(parts)
    return "|".join(sorted(tokens))


def first_token(value: str | None) -> str:
    """Return the first token from a pipe-delimited ``value``."""

    if isinstance(value, str) and value:
        return value.split("|")[0]
    return ""
