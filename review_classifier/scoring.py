"""Scoring logic for review classification."""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Tuple

from .constants import (
    MESH_EXPERIMENTAL_QUALIFIERS,
    MESH_NONREVIEW_STRONG,
    MESH_REVIEW_MARKERS,
    NONREVIEW_PUBLICATION_TYPES,
    REVIEW_PUBLICATION_TYPES,
    SRC_WEIGHTS,
)


@dataclass
class Record:
    """Normalized information required for classification."""

    ids: Mapping[str, str]
    pts: Mapping[str, List[str]]
    mesh_desc: Mapping[str, List[str]]
    mesh_qual: Mapping[str, List[str]]


@dataclass
class ScoreResult:
    """Outcome of scoring for a record."""

    label: str
    score_review: float
    score_non_review: float
    reason: str
    contrib: Dict[str, List[Tuple[str, str, str, float]]] = field(default_factory=dict)


def score_record(record: Record, logger: logging.Logger) -> ScoreResult:
    """Score a normalized record and return the classification label.

    Parameters
    ----------
    record:
        Normalized record to be classified.
    logger:
        Logger for emitting decision trail information.

    Returns
    -------
    ScoreResult
        Dataclass containing the decision and detailed scores.
    """
    contrib: Dict[str, List[Tuple[str, str, str, float]]] = {"review": [], "non_review": []}
    s_review = 0.0
    s_non = 0.0

    # Publication type contribution
    for src, pts in record.pts.items():
        weight = SRC_WEIGHTS.get(src, 0.5)
        for pt in pts:
            if pt in REVIEW_PUBLICATION_TYPES:
                val = REVIEW_PUBLICATION_TYPES[pt] * weight
                s_review += val
                contrib["review"].append(("PT", src, pt, val))
            if pt in NONREVIEW_PUBLICATION_TYPES:
                val = NONREVIEW_PUBLICATION_TYPES[pt] * weight
                s_non += val
                contrib["non_review"].append(("PT", src, pt, val))

    # Merge MeSH across sources
    desc = {d for lst in record.mesh_desc.values() for d in lst}
    qual = {q for lst in record.mesh_qual.values() for q in lst}

    for d in desc:
        if d in MESH_REVIEW_MARKERS:
            val = float(MESH_REVIEW_MARKERS[d])
            s_review += val
            contrib["review"].append(("MeSHD", "any", d, val))
        if d in MESH_NONREVIEW_STRONG:
            val = float(MESH_NONREVIEW_STRONG[d])
            s_non += val
            contrib["non_review"].append(("MeSHD", "any", d, val))

    for q in qual:
        if q in MESH_EXPERIMENTAL_QUALIFIERS:
            val = 2.0
            s_non += val
            contrib["non_review"].append(("MeSHQ", "any", q, val))

    has_meta = any("META_ANALYSIS" in pts for pts in record.pts.values())
    has_exp = len(qual.intersection(MESH_EXPERIMENTAL_QUALIFIERS)) > 0

    if has_meta and not has_exp:
        label = "review"
        reason = "PT META_ANALYSIS without experimental MeSH"
    else:
        if s_review - s_non >= 3.0:
            label = "review"
            reason = "score_gap>=3"
        elif s_non - s_review >= 3.0:
            label = "non_review"
            reason = "score_gap>=3"
        else:
            label = "uncertain"
            reason = "scores_too_close"

    result = ScoreResult(
        label=label,
        score_review=s_review,
        score_non_review=s_non,
        reason=reason,
        contrib=contrib,
    )

    logger.info(
        {
            "event": "scoring",
            "ids": record.ids,
            "score_review": s_review,
            "score_non_review": s_non,
            "label": label,
            "reason": reason,
            "contrib": sorted(
                contrib["review"] + contrib["non_review"], key=lambda x: -x[3]
            )[:10],
        }
    )

    return result
