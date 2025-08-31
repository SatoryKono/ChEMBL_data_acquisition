"""Voting logic for document classification."""
from __future__ import annotations

from typing import Dict, Iterable, Literal

from .signals import Signal
from .weights import DEFAULT_THRESHOLD

Scores = Dict[Literal["review", "nonreview"], float]


def score(signals: Iterable[Signal]) -> Scores:
    """Aggregate points from all signals into review/non-review scores."""
    totals: Scores = {"review": 0.0, "nonreview": 0.0}
    for sig in signals:
        if sig.kind == "review":
            totals["review"] += sig.points
        else:
            totals["nonreview"] += sig.points
    return totals


def decide(scores: Scores, threshold: float = DEFAULT_THRESHOLD) -> str:
    """Decide final label based on scores and threshold."""
    diff = scores["review"] - scores["nonreview"]
    if diff >= threshold:
        return "review"
    if -diff >= threshold:
        return "non-review"
    return "unknown"


__all__ = ["score", "decide"]
