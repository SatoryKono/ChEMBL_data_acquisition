"""Compatibility layer aggregating ChEMBL helpers.

This module re-exports selected public functions from :mod:`chembl_assay` and
:mod:`chembl_target` as a convenient façade.  The functions are imported
explicitly rather than via ``import *`` to make the exported API clear and to
keep linters happy.
"""

from __future__ import annotations

from ..chembl_assay import (
    get_activities,
    get_assay,
    get_assays,
    get_testitem,
)
from ..chembl_document import get_documents
from ..chembl_target import (
    extend_target,
    get_target,
    get_targets,
)
from ..clients.chembl_client import _chunked

__all__ = [
    "get_assay",
    "get_assays",
    "get_activities",
    "get_testitem",
    "get_documents",
    "get_target",
    "get_targets",
    "extend_target",
    "_chunked",
]
