"""Declarative schema package for document metadata resources."""

from __future__ import annotations

from pathlib import Path
from typing import Final

PACKAGE_DIR: Final[Path] = Path(__file__).resolve().parent
DEFAULT_DOCUMENT_SCHEMA: Final[Path] = PACKAGE_DIR / "document.yaml"

__all__ = ["DEFAULT_DOCUMENT_SCHEMA", "PACKAGE_DIR"]
