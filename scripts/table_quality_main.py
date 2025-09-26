"""Convenience wrapper matching :mod:`scripts.table_quality`."""

from __future__ import annotations

from .table_quality import main

__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
