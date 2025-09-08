"""Tests for :mod:`library.document_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library import document_postprocessing as dp


def test_postprocess_documents_creates_flags_and_sorts() -> None:
    """``postprocess_documents`` sorts rows and detects review articles."""

    path = Path("tests/data/documents_postprocess.csv")
    df = pd.read_csv(path, dtype=str)
    result = dp.postprocess_documents(df)

    # Rows sorted by computed date_code
    assert result["document_chembl_id"].tolist() == ["DOC2", "DOC1"]
    assert result["date_code"].tolist() == ["2018-08-15", "2020-05-12"]
    # Index column is zero-padded
    assert result["Index"].tolist() == ["0000", "0001"]
    # Review detection from multiple sources
    assert result["PubMed.is_review"].tolist() == [False, True]
    assert result["scholar.is_review"].tolist() == [False, True]
    assert result["OpenAlex.is_review"].tolist() == [False, True]
    # Original publication type columns have been removed
    for col in [
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
    ]:
        assert col not in result.columns
