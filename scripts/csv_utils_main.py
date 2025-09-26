"""Convenience wrapper matching :mod:`scripts.csv_utils`."""

from __future__ import annotations

from .csv_utils import main

__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - convenience entry point
    raise SystemExit(main())
