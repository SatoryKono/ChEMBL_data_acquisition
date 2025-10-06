"""Backward compatible wrapper around :mod:`library.bootstrap`."""

from __future__ import annotations

from library.bootstrap import bootstrap_cli, ensure_project_root, resolve_project_root

__all__ = ["bootstrap_cli", "ensure_project_root", "resolve_project_root"]
