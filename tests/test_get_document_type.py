"""Tests for :mod:`scripts.get_document_type` utilities."""

from __future__ import annotations

import pandas as pd

from library.utils.cli_tools.get_document_type import classify_dataframe


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

    result = classify_dataframe(df)

    assert result.loc[0, "class_label"] == "review"
