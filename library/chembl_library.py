"""Compatibility layer aggregating ChEMBL helpers."""

from __future__ import annotations

from . import chembl_assay as _chembl_assay
from . import chembl_target as _chembl_target
from .chembl_assay import *  # noqa: F401,F403
from .chembl_client import _chunked
from .chembl_target import *  # noqa: F401,F403

__all__ = [*_chembl_target.__all__, *_chembl_assay.__all__, "_chunked"]
