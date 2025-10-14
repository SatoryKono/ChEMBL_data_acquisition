"""Convenience helpers exposing canonical project paths."""

from __future__ import annotations

from pathlib import Path

__all__ = ["ROOT", "OUTPUT_DIR"]

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "data" / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
