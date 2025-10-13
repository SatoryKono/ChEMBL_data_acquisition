"""Shared filesystem paths for configuration resources."""

from __future__ import annotations

from pathlib import Path
from typing import Final

CONFIG_DIR: Final[Path] = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH: Final[Path] = CONFIG_DIR / "config.yaml"
DICTIONARY_DIR: Final[Path] = CONFIG_DIR / "dictionary"
SCHEMA_DIR: Final[Path] = CONFIG_DIR / "schema"
PIPELINE_DIR: Final[Path] = CONFIG_DIR / "pipeline"

__all__: list[str] = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DICTIONARY_DIR",
    "PIPELINE_DIR",
    "SCHEMA_DIR",
]
