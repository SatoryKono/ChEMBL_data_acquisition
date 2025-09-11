from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library import pubmed_library as pl


def test_pubmed_library_main_smoke(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run the CLI with stubbed network calls."""
    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "out.csv"

    def fake_fetch_pubmed_batch(session, pmids, delay, cfg=None):
        return [
            {
                "PubMed.PMID": pid,
                "PubMed.DOI": "10.1000/example",
                "PubMed.ArticleTitle": "Example",
            }
            for pid in pmids
        ]

    def fake_fetch_semantic_scholar_batch(session, pmids, delay, cfg=None):
        return [
            {
                "scholar.PMID": pid,
                "scholar.DOI": "10.1000/example",
                "scholar.Error": "",
            }
            for pid in pmids
        ]

    def fake_fetch_openalex(session, pmid, cfg=None, limiter=None):
        return {"OpenAlex.Error": ""}

    def fake_fetch_crossref(session, doi, cfg=None, limiter=None):
        return {"crossref.Error": ""}

    class DummyLimiter:
        def acquire(self) -> None:
            return None

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        pl, "fetch_semantic_scholar_batch", fake_fetch_semantic_scholar_batch
    )
    monkeypatch.setattr(pl, "fetch_openalex", fake_fetch_openalex)
    monkeypatch.setattr(pl, "fetch_crossref", fake_fetch_crossref)
    monkeypatch.setattr(pl, "get_limiter", lambda *args, **kwargs: DummyLimiter())

    exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty
