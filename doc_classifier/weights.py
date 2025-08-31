"""Weight configuration for document classification."""
from __future__ import annotations

from typing import Dict

Weights = Dict[str, float]
ScoreSchema = Dict[str, Dict[str, float]]

SOURCE_WEIGHTS: Weights = {
    "pubmed": 3.0,
    "openalex": 2.0,
    "crossref": 1.5,
    "scholar": 1.0,
}

SCORE_SCHEMA: ScoreSchema = {
    "PT": {"review": 2.0, "non-review": 1.5},
    "MeSH_desc": {"review": 1.0, "non-review": 0.8},
    "MeSH_qual": {"review": 0.5, "non-review": 0.6},
}

DEFAULT_THRESHOLD: float = 1.0

__all__ = ["SOURCE_WEIGHTS", "SCORE_SCHEMA", "DEFAULT_THRESHOLD"]
