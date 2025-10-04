"""Shared filesystem paths for configuration resources."""
from __future__ import annotations

from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = CONFIG_DIR / "config.yaml"
DICTIONARY_DIR = CONFIG_DIR / "dictionary"

__all__ = ["CONFIG_DIR", "DEFAULT_CONFIG_PATH", "DICTIONARY_DIR"]
