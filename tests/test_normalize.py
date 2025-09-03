
from doc_classifier.terms import parse_terms

def test_split_and_normalize():
    assert split_and_normalize("Review|Journal Article", "|;,/") == ["review", "journal article"]


def test_split_and_canon_basic():
    text = "Review Article; Meta-Analysis / randomized Controlled Trial"
    tokens = parse_terms(text)
    assert tokens == [
        "meta-analysis",
        "randomized controlled trial",
        "review",
    ]

def test_is_review():
    assert is_review({"review"})
    assert is_review({"systematic review"})
    assert not is_review({"journal article", "letter"})


def test_split_and_canon_synonyms():
    text = "Mini Review | SLR, JournalArticle"
    tokens = parse_terms(text)
    assert tokens == ["journal-article", "review", "systematic review"]
