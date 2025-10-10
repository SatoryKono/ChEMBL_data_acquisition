"""Configuration package exposing shared file-system paths."""

from __future__ import annotations

from .paths import (
    CONFIG_DIR,
    DEFAULT_CONFIG_PATH,
    DICTIONARY_DIR,
    PIPELINE_DIR,
    SCHEMA_DIR,
)

__all__ = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DICTIONARY_DIR",
    "PIPELINE_DIR",
    "SCHEMA_DIR",
]
