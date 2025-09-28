"""Convenience wrapper matching :mod:`scripts.mapper`."""

from __future__ import annotations

from .mapper import main

__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
