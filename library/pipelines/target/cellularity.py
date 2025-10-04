"""Wrappers around organism cellularity helpers used by the targets pipeline."""

from __future__ import annotations

from .organism_classification import (
    DEFAULT_RULES,
    OrganismClassificationRules,
    add_cellularity,
    add_cellularity_smart,
    classify_by_lineage,
    classify_record,
    normalize,
)

__all__ = [
    "DEFAULT_RULES",
    "OrganismClassificationRules",
    "add_cellularity",
    "add_cellularity_smart",
    "classify_by_lineage",
    "classify_record",
    "normalize",
]
