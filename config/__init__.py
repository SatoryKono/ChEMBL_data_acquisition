"""Configuration package exposing shared file-system paths."""
from __future__ import annotations

from .paths import (
    CONFIG_DIR,
    CONFIG_SCHEMA_PATH,
    DEFAULT_CONFIG_PATH,
    DICTIONARY_DIR,
    POSTPROCESSING_CONFIG_DIR,
    SCHEMA_DIR,
)

__all__ = [
    "CONFIG_DIR",
    "CONFIG_SCHEMA_PATH",
    "DEFAULT_CONFIG_PATH",
    "DICTIONARY_DIR",
    "POSTPROCESSING_CONFIG_DIR",
    "SCHEMA_DIR",
]
