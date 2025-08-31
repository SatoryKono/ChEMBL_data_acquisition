"""Constants and mappings for review classification."""
from __future__ import annotations

from typing import Dict, List, Mapping

# Source weights compensate for varying accuracy of publication type fields
SRC_WEIGHTS: Mapping[str, float] = {
    "pubmed": 1.0,
    "crossref": 0.7,
    "openalex": 0.7,
    "scholar": 0.5,
}

# Publication types signalling a review
REVIEW_PUBLICATION_TYPES: Mapping[str, int] = {
    "REVIEW": 3,
    "SYSTEMATIC_REVIEW": 4,
    "META_ANALYSIS": 5,
    "SCOPING_REVIEW": 3,
}

# Publication types indicating non-review content
NONREVIEW_PUBLICATION_TYPES: Mapping[str, int] = {
    "RCT": 3,
    "CLINICAL_TRIAL": 3,
    "CASE_REPORT": 3,
    "PRECLINICAL": 3,
    "IN_VITRO": 3,
    "IN_VIVO": 3,
}

# MeSH descriptor markers for review-like content
MESH_REVIEW_MARKERS: Mapping[str, int] = {
    "review literature as topic": 4,
    "systematic reviews as topic": 4,
    "meta-analysis as topic": 4,
    "practice guidelines as topic": 4,
    "evidence-based practice": 4,
}

# Strong MeSH descriptors indicating experimental/non-review work
MESH_NONREVIEW_STRONG: Mapping[str, int] = {
    "in vitro techniques": 3,
    "cell line": 3,
    "animals": 3,
    "disease models, animal": 3,
    "randomized controlled trial as topic": 3,
    "double-blind method": 3,
    "treatment outcome": 3,
}

# MeSH qualifiers signalling experimental work
MESH_EXPERIMENTAL_QUALIFIERS: List[str] = [
    "methods",
    "drug therapy",
    "adverse effects",
    "chemistry",
    "metabolism",
    "pharmacology",
    "radiation effects",
    "ultrastructure",
    "genetics",
    "immunology",
]

# Mapping of source specific publication type labels to normalized ones
PUBLICATION_TYPE_MAPPING: Dict[str, Dict[str, str]] = {
    "pubmed": {
        "Review": "REVIEW",
        "Systematic Review": "SYSTEMATIC_REVIEW",
        "Meta-Analysis": "META_ANALYSIS",
        "Scoping Review": "SCOPING_REVIEW",
        "Clinical Trial": "CLINICAL_TRIAL",
        "Randomized Controlled Trial": "RCT",
        "Case Reports": "CASE_REPORT",
    },
    "crossref": {
        "review": "REVIEW",
        "meta-analysis": "META_ANALYSIS",
        "journal-article": "OTHER",
    },
    "openalex": {
        "review": "REVIEW",
        "meta_analysis": "META_ANALYSIS",
        "clinical_trial": "CLINICAL_TRIAL",
    },
    "scholar": {
        "Review": "REVIEW",
        "Meta Analysis": "META_ANALYSIS",
    },
}
