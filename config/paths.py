"""Shared filesystem paths for configuration resources."""
from __future__ import annotations

from pathlib import Path

CONFIG_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = CONFIG_DIR / "config.yaml"
DICTIONARY_DIR = CONFIG_DIR / "dictionary"
SCHEMA_DIR = CONFIG_DIR / "schema"
POSTPROCESSING_CONFIG_DIR = CONFIG_DIR / "postprocessing"
CONFIG_SCHEMA_PATH = CONFIG_DIR.parent / "config.schema.json"

__all__ = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DICTIONARY_DIR",
    "SCHEMA_DIR",
    "POSTPROCESSING_CONFIG_DIR",
    "CONFIG_SCHEMA_PATH",
]
