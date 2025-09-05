"""Tests for :mod:`library.pubmed_library`."""

from __future__ import annotations

from pathlib import Path

from library import pubmed_library as pl


DATA_DIR = Path(__file__).parent / "data"


def test_read_pmids() -> None:
    path = DATA_DIR / "pmids.csv"
    pmids = pl.read_pmids(path)
    assert pmids == ["1", "2"]
