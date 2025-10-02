"""Tests for :mod:`library.document_postprocessing`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library import document_postprocessing as dp
from library.config import IoCfg


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


def test_postprocess_file_roundtrip(tmp_path: Path) -> None:
    """``postprocess_file`` writes metadata and normalises NA handling."""

    df = pd.DataFrame(
        {
            "document_chembl_id": ["DOC1", "DOC2"],
            "title": ["Example", None],
            "PubMed.PublicationType": ["Review", ""],
            "scholar.PublicationTypes": ["Review", None],
            "OpenAlex.PublicationTypes": ["review-article", ""],
            "Index": [0, 1],
        }
    )
    cfg = IoCfg(csv_sep=";", csv_encoding="utf8")
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False, sep=cfg.csv_sep, encoding=cfg.csv_encoding)
    output_path = tmp_path / "out.csv"

    dp.postprocess_file(input_path, output_path, cfg=cfg)

    result = pd.read_csv(
        output_path,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        keep_default_na=False,
    )
    assert result["title"].tolist() == ["Example", "nan"]
    assert not result.isna().any().any()
    assert result["PubMed.MonthCompleted"].tolist() == ["", ""]

    meta_path = Path(f"{output_path}.meta.yaml")
    assert meta_path.exists()
