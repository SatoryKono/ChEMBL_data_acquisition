"""Unit tests for review classification pipeline."""
from __future__ import annotations

import logging

from review_classifier.normalization import normalize_publication_types
from review_classifier.scoring import Record, score_record


def test_normalize_publication_types() -> None:
    raw = {"pubmed": ["Systematic Review"], "crossref": ["meta-analysis", "journal-article"]}
    norm = normalize_publication_types(raw)
    assert norm["pubmed"] == ["SYSTEMATIC_REVIEW"]
    assert "META_ANALYSIS" in norm["crossref"]


def test_score_record_review() -> None:
    logger = logging.getLogger("test")
    logger.addHandler(logging.NullHandler())
    record = Record(
        ids={"pmid": "1"},
        pts={"pubmed": ["SYSTEMATIC_REVIEW"]},
        mesh_desc={"pubmed": ["systematic reviews as topic"]},
        mesh_qual={},
    )
    result = score_record(record, logger)
    assert result.label == "review"


def test_score_record_non_review() -> None:
    logger = logging.getLogger("test")
    logger.addHandler(logging.NullHandler())
    record = Record(
        ids={"pmid": "2"},
        pts={"pubmed": ["CLINICAL_TRIAL"], "openalex": ["REVIEW"]},
        mesh_desc={"pubmed": ["double-blind method"]},
        mesh_qual={"pubmed": ["methods"]},
    )
    result = score_record(record, logger)
    assert result.label == "non_review"
