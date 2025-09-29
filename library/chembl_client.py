"""Backward-compatible wrapper for :mod:`library.clients.chembl`."""

from __future__ import annotations

from .clients.chembl import ChemblClient, _chunked

__all__ = ["ChemblClient", "_chunked"]
