from doc_classifier.classifier import compute_scores, decide_label
from doc_classifier.terms import parse_terms


def classify(pubmed: str, scholar: str, openalex: str):
    scores = compute_scores(
        parse_terms(pubmed), parse_terms(scholar), parse_terms(openalex)
    )
    label = decide_label(scores)
    return label, scores


def test_review_from_pubmed():
    label, scores = classify("Journal Article|Review", "", "")
    assert label == "review"
    assert scores["review"] > 0


def test_unknown_editorial():
    label, scores = classify("Editorial", "", "")
    assert label == "unknown"
    assert scores["unknown"] > 0


def test_experimental_proceedings():
    label, scores = classify("", "", "proceedings-article")
    assert label == "experimental"
    assert scores["experimental"] > 0


def test_unknown_empty():
    label, scores = classify("", "", "")
    assert label == "unknown"
    assert scores == {"review": 0, "experimental": 0, "unknown": 0}


def test_conflict_review_and_posted():
    label, scores = classify("Review", "", "posted-content")
    assert label == "review"
