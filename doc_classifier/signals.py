"""Extraction of review/non-review signals from normalised tokens."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Literal

from .weights import SCORE_SCHEMA, SOURCE_WEIGHTS


@dataclass
class Signal:
    """Single classification signal fired from an input field."""

    source: Literal["pubmed", "openalex", "crossref", "scholar"]
    field: Literal["PT", "MeSH_desc", "MeSH_qual"]
    kind: Literal["review", "non-review"]
    token: str
    weight: float
    points: float


# ---------------------------------------------------------------------------
# Signal dictionaries
# ---------------------------------------------------------------------------

# PublicationType signals
REVIEW_PT = {
    "review",
    "systematic review",
    "meta-analysis",
    "scoping review",
    "umbrella review",
}

NONREVIEW_PT = {
    "clinical trial",
    "randomized controlled trial",
    "controlled clinical trial",
    "case report",
    "case series",
    "cohort study",
    "case-control study",
    "cross-sectional study",
    "longitudinal study",
    "comparative study",
    "multicenter study",
    "evaluation study",
    "validation study",
    "replication study",
    "observational study",
    "letter",
    "short report",
    "technical note",
    "methods paper",
}

# MeSH descriptor signals
REVIEW_DESC = {
    "systematic reviews as topic",
    "meta-analysis as topic",
    "review",
    "evidence-based practice",
    "evidence synthesis",
}

NONREVIEW_DESC = {
    "drug evaluation",
    "preclinical",
    "in vitro techniques",
    "cell line",
    "disease models",
    "animal",
    "animals",
    "mice",
    "rats",
    "humans",
    "tissue culture techniques",
    "flow cytometry",
    "spectrometry",
    "mass",
    "chromatography",
    "high pressure liquid",
    "x-ray crystallography",
    "western blotting",
    "microscopy",
    "electron",
}

# MeSH qualifier signals
REVIEW_QUAL = {"review", "analysis", "methods"}

NONREVIEW_QUAL = {
    "therapy",
    "drug effects",
    "metabolism",
    "chemistry",
    "pharmacology",
    "physiology",
    "pathology",
    "enzymology",
    "antagonists & inhibitors",
    "administration & dosage",
    "adverse effects",
    "diagnosis",
    "genetics",
    "immunology",
}


# ---------------------------------------------------------------------------
# Extraction helpers
# ---------------------------------------------------------------------------


def _build_signal(source: str, field: str, kind: str, token: str) -> Signal:
    weight = SOURCE_WEIGHTS[source]
    points = weight * SCORE_SCHEMA[field][kind]
    return Signal(source=source, field=field, kind=kind, token=token, weight=weight, points=points)


def extract_pt_signals(tokens: List[str], source: str) -> List[Signal]:
    """Extract signals from PublicationType tokens."""
    signals: List[Signal] = []
    for token in tokens:
        if token in REVIEW_PT:
            signals.append(_build_signal(source, "PT", "review", token))
        elif token in NONREVIEW_PT:
            signals.append(_build_signal(source, "PT", "non-review", token))
    return signals


def extract_mesh_signals(
    desc: List[str], qual: List[str], source: str
) -> List[Signal]:
    """Extract signals from MeSH descriptor and qualifier tokens."""
    signals: List[Signal] = []
    for token in desc:
        if token in REVIEW_DESC:
            signals.append(_build_signal(source, "MeSH_desc", "review", token))
        elif token in NONREVIEW_DESC:
            signals.append(_build_signal(source, "MeSH_desc", "non-review", token))
    for token in qual:
        if token in REVIEW_QUAL:
            signals.append(_build_signal(source, "MeSH_qual", "review", token))
        elif token in NONREVIEW_QUAL:
            signals.append(_build_signal(source, "MeSH_qual", "non-review", token))
    return signals


__all__ = [
    "Signal",
    "extract_pt_signals",
    "extract_mesh_signals",
    "REVIEW_PT",
    "NONREVIEW_PT",
    "REVIEW_DESC",
    "NONREVIEW_DESC",
    "REVIEW_QUAL",
    "NONREVIEW_QUAL",
]
