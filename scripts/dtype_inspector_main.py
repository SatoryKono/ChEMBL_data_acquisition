"""Convenience wrapper matching :mod:`scripts.dtype_inspector`."""

from __future__ import annotations

from .dtype_inspector import main

__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
