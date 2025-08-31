"""Transformations and classification logic for document types."""
from __future__ import annotations

from dataclasses import dataclass
import json
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

# Mapping of raw publication type strings to normalized categories
PT_MAPPING: Dict[str, str] = {
    "review": "REVIEW",
    "systematic review": "SYSTEMATIC_REVIEW",
    "meta-analysis": "META_ANALYSIS",
    "meta analysis": "META_ANALYSIS",
    "scoping review": "SCOPING_REVIEW",
    "guideline": "GUIDELINE",
    "randomized controlled trial": "RCT",
    "randomised controlled trial": "RCT",
    "clinical trial": "CLINICAL_TRIAL",
    "case report": "CASE_REPORT",
    "case reports": "CASE_REPORT",
    "preclinical": "PRECLINICAL",
    "in vitro": "IN_VITRO",
    "in vivo": "IN_VIVO",
    "protocol": "PROTOCOL",
    "editorial": "EDITORIAL",
    "letter": "LETTER",
    "data paper": "DATA_PAPER",
}

# Score weights
SOURCE_WEIGHTS = {
    "pubmed": 1.0,
    "crossref": 0.7,
    "openalex": 0.7,
    "scholar": 0.5,
}

REVIEW_PTS = {"REVIEW", "SYSTEMATIC_REVIEW", "META_ANALYSIS", "SCOPING_REVIEW"}
NON_REVIEW_PTS = {
    "RCT",
    "CLINICAL_TRIAL",
    "CASE_REPORT",
    "PRECLINICAL",
    "IN_VITRO",
    "IN_VIVO",
}

PT_SCORES_REVIEW = {
    "REVIEW": 3,
    "SYSTEMATIC_REVIEW": 4,
    "META_ANALYSIS": 5,
    "SCOPING_REVIEW": 3,
}

PT_SCORES_NON_REVIEW = {
    "RCT": 3,
    "CLINICAL_TRIAL": 3,
    "CASE_REPORT": 3,
    "PRECLINICAL": 3,
    "IN_VITRO": 3,
    "IN_VIVO": 3,
}

REVIEW_MESH = {
    "review literature as topic",
    "systematic reviews as topic",
    "meta-analysis as topic",
    "practice guidelines as topic",
    "evidence-based practice",
}

EXPERIMENTAL_MESH_DESCRIPTORS = {
    "in vitro techniques",
    "cell line",
    "cell lines",
    "animals",
    "disease models, animal",
    "randomized controlled trial as topic",
    "double-blind method",
    "treatment outcome",
}

EXPERIMENTAL_MESH_QUALIFIERS = {
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
}

RE_RULES_VERSION = "1.0.0"


@dataclass
class ClassificationResult:
    """Result of document type classification."""

    record_id: str
    label: str
    reason: str
    score_review: float
    score_non_review: float
    factors: List[str]
    normalized_pt: Dict[str, List[str]]
    normalized_mesh: Dict[str, List[str]]


def _split(cell: str) -> List[str]:
    """Split a cell string by common delimiters and normalise case."""

    if not cell:
        return []
    parts = re.split(r"[;|,]", cell)
    return [p.strip().lower() for p in parts if p.strip()]


def _map_pt(values: Iterable[str]) -> List[str]:
    mapped = []
    for val in values:
        mapped.append(PT_MAPPING.get(val, val.upper()))
    return mapped


def classify_row(row: pd.Series, log_file) -> ClassificationResult:
    """Classify a single record.

    Parameters
    ----------
    row:
        Row of the input DataFrame.
    log_file:
        Open file handle to append NDJSON logs.

    Returns
    -------
    ClassificationResult
        Classification outcome and auxiliary information.
    """

    record_id = row.get("record_id", "")

    original_pt = {
        "pubmed": row.get("pubmed_publication_types", ""),
        "crossref": row.get("crossref_type", ""),
        "openalex": row.get("openalex_type", ""),
        "scholar": row.get("scholar_type", ""),
    }
    normalized_pt = {src: _map_pt(_split(val)) for src, val in original_pt.items()}

    original_mesh = {
        "pubmed_descriptors": row.get("pubmed_mesh_descriptors", ""),
        "pubmed_qualifiers": row.get("pubmed_mesh_qualifiers", ""),
        "openalex_descriptors": row.get("openalex_mesh_descriptors", ""),
        "openalex_qualifiers": row.get("openalex_mesh_qualifiers", ""),
    }
    normalized_mesh = {src: _split(val) for src, val in original_mesh.items()}

    # Stage 1: direct decision if review PT without conflicts
    all_pts = [pt for pts in normalized_pt.values() for pt in pts]
    if any(pt in REVIEW_PTS for pt in all_pts) and not any(
        pt in NON_REVIEW_PTS for pt in all_pts
    ):
        reason = "direct_review_pt"
        result = ClassificationResult(
            record_id,
            "review",
            reason,
            score_review=0,
            score_non_review=0,
            factors=["publication_type"],
            normalized_pt=normalized_pt,
            normalized_mesh=normalized_mesh,
        )
        log_file.write(json.dumps({"record_id": record_id, "label": "review", "reason": reason, "pt": original_pt, "mesh": original_mesh}) + "\n")
        return result

    score_review = 0.0
    score_non_review = 0.0
    factors: List[str] = []

    # Publication type scoring
    for src, pts in normalized_pt.items():
        weight = SOURCE_WEIGHTS.get(src, 1.0)
        for pt in pts:
            if pt in PT_SCORES_REVIEW:
                inc = PT_SCORES_REVIEW[pt] * weight
                score_review += inc
                factors.append(f"{src}:{pt}:{inc}")
            if pt in PT_SCORES_NON_REVIEW:
                inc = PT_SCORES_NON_REVIEW[pt] * weight
                score_non_review += inc
                factors.append(f"{src}:{pt}:-{inc}")

    # MeSH scoring (no source weights specified)
    descriptors = normalized_mesh["pubmed_descriptors"] + normalized_mesh["openalex_descriptors"]
    qualifiers = normalized_mesh["pubmed_qualifiers"] + normalized_mesh["openalex_qualifiers"]

    for d in descriptors:
        if d in REVIEW_MESH:
            score_review += 4
            factors.append(f"mesh_descriptor:{d}:4")
        if d in EXPERIMENTAL_MESH_DESCRIPTORS:
            score_non_review += 3
            factors.append(f"mesh_descriptor:{d}:-3")
    for q in qualifiers:
        if q in EXPERIMENTAL_MESH_QUALIFIERS:
            score_non_review += 2
            factors.append(f"mesh_qualifier:{q}:-2")

    # Special rules
    if "META_ANALYSIS" in all_pts and not any(q in EXPERIMENTAL_MESH_QUALIFIERS for q in qualifiers):
        reason = "meta_analysis_no_experiment"
        label = "review"
    elif "GUIDELINE" in all_pts and "practice guidelines as topic" in descriptors and not any(
        d in EXPERIMENTAL_MESH_DESCRIPTORS for d in descriptors
    ):
        reason = "guideline_with_topic"
        label = "review"
    elif "PROTOCOL" in all_pts and (
        "systematic review protocol" in row.get("title", "").lower()
        or any(d in REVIEW_MESH for d in descriptors)
    ):
        reason = "protocol_review"
        label = "review"
    else:
        if score_review - score_non_review >= 3:
            label = "review"
            reason = "score_threshold"
        elif score_non_review - score_review >= 3:
            label = "non_review"
            reason = "score_threshold"
        else:
            label = "uncertain"
            reason = "low_confidence"

    if not factors:
        label = "uncertain"
        reason = "no_signals"

    top_factors = sorted(factors, key=lambda x: float(x.split(":")[-1]), reverse=True)[:5]

    result = ClassificationResult(
        record_id,
        label,
        reason,
        score_review,
        score_non_review,
        top_factors,
        normalized_pt,
        normalized_mesh,
    )

    log_file.write(
        json.dumps(
            {
                "record_id": record_id,
                "original_pt": original_pt,
                "original_mesh": original_mesh,
                "normalized_pt": normalized_pt,
                "normalized_mesh": normalized_mesh,
                "score_review": score_review,
                "score_non_review": score_non_review,
                "factors": factors,
                "label": label,
                "reason": reason,
            }
        )
        + "\n"
    )

    return result
