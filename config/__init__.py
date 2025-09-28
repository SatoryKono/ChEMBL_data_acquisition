"""Configuration package exposing shared file-system paths."""

from __future__ import annotations

from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = CONFIG_DIR / "config.yaml"

__all__ = ["CONFIG_DIR", "DEFAULT_CONFIG_PATH"]
