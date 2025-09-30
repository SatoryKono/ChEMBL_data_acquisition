from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library import pubmed_library as pl
from library.pubmed import aggregation as pa


def test_pubmed_library_main_smoke(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Run the CLI with stubbed network calls."""
    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "out.csv"
    verbose_output_csv = tmp_path / "out_verbose.csv"

    def fake_fetch_pubmed_batch(session, pmids, delay, cfg=None, *, client=None):
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

    levels: list[str] = []
    dumps: list[str] = []

    def fake_emit_output(log, level: str, output: str) -> None:
        levels.append(level)
        dumps.append(output)

    monkeypatch.setattr(pa, "_emit_output", fake_emit_output)

    exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "INFO",
        ]
    )
    assert exit_code == 0
    assert levels
    default_levels = levels.copy()
    default_dumps = dumps.copy()
    assert all(level == "DEBUG" for level in default_levels)
    assert default_dumps
    assert all(len(dump) < 200 for dump in default_dumps)
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty

    start_idx = len(levels)
    start_dump_idx = len(dumps)
    exit_code_verbose = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(verbose_output_csv),
            "--log-level",
            "INFO",
            "--verbose",
        ]
    )
    assert exit_code_verbose == 0
    verbose_levels = levels[start_idx:]
    verbose_dumps = dumps[start_dump_idx:]
    assert verbose_levels
    assert all(level == "INFO" for level in verbose_levels)
    assert verbose_dumps
    for dump in verbose_dumps:
        assert "PMID" in dump
        assert "DOI" in dump
        assert "Title" in dump
        assert "PubMed.PMID" not in dump
        assert len(dump) < 200
