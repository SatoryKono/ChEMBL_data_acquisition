"""Thin wrappers around CrossRef APIs."""

from .crossref import fetch_crossref
from .crossref_client import CrossrefClient

__all__ = ["fetch_crossref", "CrossrefClient"]
