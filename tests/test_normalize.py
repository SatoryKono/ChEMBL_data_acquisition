import pytest

from doc_classifier.normalize import split_and_canon


def test_split_and_canon_basic():
    text = "Review Article; Meta-Analysis / randomized Controlled Trial"
    tokens = split_and_canon(text)
    assert tokens == [
        "review",
        "meta-analysis",
        "randomized controlled trial",
    ]


def test_split_and_canon_synonyms():
    text = "Mini Review | SLR, JournalArticle"
    tokens = split_and_canon(text)
    assert tokens == ["review", "systematic review", "journal-article"]
