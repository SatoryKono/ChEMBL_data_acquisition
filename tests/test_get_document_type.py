"""Tests for :mod:`get_document_type` utilities."""

from __future__ import annotations

import pandas as pd

from get_document_type import classify_dataframe


def test_classify_dataframe_basic() -> None:
    """Basic smoke test for :func:`classify_dataframe`."""

    df = pd.DataFrame(
        {
            "title": ["t"],
            "abstract": ["a"],
            "PubMed.PublicationType": ["review"],
            "scholar.PublicationTypes": [""],
            "OpenAlex.PublicationTypes": [""],
            "OpenAlex.TypeCrossref": [""],
        }
    )

    result = classify_dataframe(
        df,
        min_review_score=1,
        min_unknown_score=1,
        min_experimental_score=1,
    )

    assert result.loc[0, "class_label"] == "review"
