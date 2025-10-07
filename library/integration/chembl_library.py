"""Compatibility layer aggregating ChEMBL helpers.

This module re-exports selected public functions from :mod:`chembl_assay` and
:mod:`chembl_target` as a convenient façade.  The functions are imported
explicitly rather than via ``import *`` to make the exported API clear and to
keep linters happy.
"""

from __future__ import annotations

from library.clients import _chunked

from ..pipelines.assay import get_activities, get_assay, get_assays, get_testitem
from ..pipelines.document import get_documents
from ..pipelines.target.chembl_target import (
    extend_target,
    get_target,
    get_target_payload,
    get_targets,
    get_targets_payloads,
    get_targets_raw_frame,
    iter_target_batches,
    iter_target_batches_with_retry,
)
from ..pipelines.tissue import get_tissues

__all__ = [
    "get_assay",
    "get_assays",
    "get_activities",
    "get_tissues",
    "get_testitem",
    "get_documents",
    "get_target",
    "get_target_payload",
    "get_targets",
    "get_targets_payloads",
    "get_targets_raw_frame",
    "iter_target_batches",
    "iter_target_batches_with_retry",
    "extend_target",
    "_chunked",
]
