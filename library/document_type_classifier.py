"""Scoring and decision logic for publication classification."""
from __future__ import annotations

from typing import Dict, Iterable, Mapping

from .document_type_terms import REVIEW_TERMS, EXPERIMENTAL_TERMS, UNKNOWN_TERMS

# Source weights used in weighted voting
SOURCE_WEIGHTS: Mapping[str, int] = {
    "pubmed": 4,
    "openalex": 3,
    "scholar": 2,
}


def compute_scores(
    pubmed_terms: Iterable[str],
    scholar_terms: Iterable[str],
    openalex_terms: Iterable[str],
) -> Dict[str, int]:
    """Calculate scores for each class based on detected terms.

    Parameters
    ----------
    pubmed_terms, scholar_terms, openalex_terms:
        Normalised term lists from the respective sources.

    Returns
    -------
    dict[str, int]
        Mapping of class names (``review``, ``experimental``, ``unknown``) to
        accumulated weights.
    """
    scores = {"review": 0, "experimental": 0, "unknown": 0}

    def add_scores(terms: Iterable[str], weight: int) -> None:
        for term in terms:
            if term in REVIEW_TERMS:
                scores["review"] += weight
            if term in EXPERIMENTAL_TERMS:
                scores["experimental"] += weight
            if term in UNKNOWN_TERMS:
                scores["unknown"] += weight

    add_scores(pubmed_terms, SOURCE_WEIGHTS["pubmed"])
    add_scores(scholar_terms, SOURCE_WEIGHTS["scholar"])
    add_scores(openalex_terms, SOURCE_WEIGHTS["openalex"])
    return scores


def decide_label(
    scores: Mapping[str, int],
    *,
    min_review_score: int = 1,
    min_unknown_score: int = 2,
    min_experimental_score: int = 1,
) -> str:
    """Decide final class label based on score comparison.

    Parameters
    ----------
    scores:
        Mapping with keys ``review``, ``experimental`` and ``unknown``.
    min_review_score, min_unknown_score, min_experimental_score:
        Minimum scores required for the respective classes.

    Returns
    -------
    str
        One of ``review``, ``experimental`` or ``unknown``.
    """
    r = scores.get("review", 0)
    e = scores.get("experimental", 0)
    u = scores.get("unknown", 0)

    if r > max(e, u) and r >= min_review_score:
        return "review"
    if u > max(r, e) and u >= min_unknown_score:
        return "unknown"
    if e > max(r, u) and e >= min_experimental_score:
        return "experimental"
    if r == 0 and e == 0 and u > 0:
        return "unknown"
    return "unknown"
