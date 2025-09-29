from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import library.clients.pubmed as pl


class DummyLimiter:
    def acquire(self) -> None:  # pragma: no cover - simple stub
        return None


def _stub_network(monkeypatch: pytest.MonkeyPatch) -> None:
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

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        pl, "fetch_semantic_scholar_batch", fake_fetch_semantic_scholar_batch
    )
    monkeypatch.setattr(pl, "fetch_openalex", lambda *a, **k: {"OpenAlex.Error": ""})
    monkeypatch.setattr(pl, "fetch_crossref", lambda *a, **k: {"crossref.Error": ""})
    monkeypatch.setattr(pl, "get_limiter", lambda *a, **k: DummyLimiter())


def test_pubmed_single_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure :func:`library.clients.pubmed.main` creates ``Config`` only once."""
    _stub_network(monkeypatch)
    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "out.csv"

    calls: dict[str, int] = {"cfg": 0}
    OriginalConfig = pl.Config

    class CountingConfig(OriginalConfig):
        def __init__(self, *args, **kwargs):
            calls["cfg"] += 1
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(pl, "Config", CountingConfig)

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
    assert calls["cfg"] == 1
