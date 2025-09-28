"""Utilities shared across API client implementations."""

from __future__ import annotations

from collections.abc import Iterable, Iterator


def chunked(items: Iterable[str], size: int) -> Iterator[list[str]]:
  """Yield ``size``-sized lists from ``items``."""

  if size <= 0:
    raise ValueError("size must be a positive integer")

  chunk: list[str] = []
  for item in items:
    chunk.append(item)
    if len(chunk) == size:
      yield chunk
      chunk = []
  if chunk:
    yield chunk


__all__ = ["chunked"]
