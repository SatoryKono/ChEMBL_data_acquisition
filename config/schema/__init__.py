"""Schema package exposing packaged YAML resources."""

from __future__ import annotations

from importlib.resources import files
from typing import Final

DOCUMENT_SCHEMA_RESOURCE = files(__name__).joinpath("document.yaml")

__all__ = ["DOCUMENT_SCHEMA_RESOURCE"]
