"""Tests for :mod:`library.document_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pytest

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


@pytest.fixture
def documents_with_numeric_date_parts() -> pd.DataFrame:
    """Return rows with numeric date components that require padding."""

    return pd.DataFrame(
        {
            "document_chembl_id": ["DOC_A", "DOC_B", "DOC_C"],
            "PubMed.DayCompleted": [2, 0, 0],
            "PubMed.MonthCompleted": [5, 9, 3],
            "PubMed.YearCompleted": [2020, 2019, 85],
            "PubMed.DayRevised": [0, 7, 0],
            "PubMed.MonthRevised": [0, 9, 0],
            "PubMed.YearRevised": [0, 2019, 0],
        }
    )


def test_postprocess_documents_zero_pads_date_components(
    documents_with_numeric_date_parts: pd.DataFrame,
) -> None:
    """Numeric month/day values are padded in the generated ``date_code``."""

    result = dp.postprocess_documents(documents_with_numeric_date_parts)
    date_by_doc = dict(zip(result["document_chembl_id"], result["date_code"]))

    assert date_by_doc["DOC_A"] == "2020-05-02"
    assert date_by_doc["DOC_B"] == "2019-09-07"
    assert date_by_doc["DOC_C"] == "0085-03-01"
