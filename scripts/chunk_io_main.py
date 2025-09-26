"""Convenience wrapper matching :mod:`scripts.chunk_io`."""

from __future__ import annotations

from .chunk_io import main

__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
