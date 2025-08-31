from doc_classifier import normalize, signals, vote


def classify_pt(value: str, source: str):
    tokens = normalize.split_and_canon(value)
    sigs = signals.extract_pt_signals(tokens, source)
    scores = vote.score(sigs)
    label = vote.decide(scores)
    return label, scores


def test_review_from_scholar():
    label, scores = classify_pt("Narrative Review", "scholar")
    assert label == "review"
    assert scores["review"] > scores["nonreview"]


def test_nonreview_from_pubmed():
    label, scores = classify_pt("Clinical Trial", "pubmed")
    assert label == "non-review"
    assert scores["nonreview"] > scores["review"]


def test_unknown_no_signal():
    label, _ = classify_pt("Journal Article", "pubmed")
    assert label == "unknown"
