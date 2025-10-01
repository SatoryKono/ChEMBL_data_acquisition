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
    verbose_output_csv = tmp_path / "out_verbose.csv"

    def fake_fetch_pubmed_batch(session, pmids, delay, cfg=None, *, client=None):
        records = []
        for index, pid in enumerate(pmids):
            record = {
                "PubMed.PMID": pid,
                "PubMed.DOI": "10.1000/example",
                "PubMed.ArticleTitle": "Example",
            }
            if index == 0:
                record.pop("PubMed.PMID")
            records.append(record)
        return records

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

    def fake_print_results(records, *, level: str = "DEBUG") -> None:
        levels.extend([level] * len(records))

    monkeypatch.setattr(pl, "print_results", fake_print_results)

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
    assert all(level == "DEBUG" for level in levels)
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty
    meta_path = Path(f"{output_csv}.meta.yaml")
    assert meta_path.exists()
    failure_path = output_csv.with_name(f"{output_csv.stem}_failure_cases.csv")
    assert failure_path.exists()
    quality_path = Path(f"{output_csv.with_suffix('')}_quality_report_table.csv")
    correlation_path = Path(
        f"{output_csv.with_suffix('')}_data_correlation_report_table.csv"
    )
    assert quality_path.exists()
    assert correlation_path.exists()
    failure_df = pd.read_csv(failure_path)
    assert not failure_df.empty

    start_idx = len(levels)
    exit_code_verbose = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(verbose_output_csv),
            "--log-level",
            "INFO",
            "--keep-verbose-dumps",
        ]
    )
    assert exit_code_verbose == 0
    verbose_levels = levels[start_idx:]
    assert verbose_levels
    assert all(level == "INFO" for level in verbose_levels)
