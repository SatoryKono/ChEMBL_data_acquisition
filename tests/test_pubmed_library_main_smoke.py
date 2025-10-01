from __future__ import annotations

from pathlib import Path
import hashlib

import pandas as pd
import pytest
import yaml

from library import pubmed_library as pl


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
    failure_path = output_csv.with_name(f"{output_csv.stem}_failure_cases.csv")
    assert not failure_path.exists()

    meta_path = output_csv.with_name(output_csv.name + ".meta.yaml")
    assert meta_path.exists()
    digest = hashlib.sha256(output_csv.read_bytes()).hexdigest()
    meta = yaml.safe_load(meta_path.read_text())
    assert meta["stats"]["output_sha256"] == digest
    assert meta["schema"] == "PubMedDocumentsSchema"

    start_idx = len(levels)
    second_output = tmp_path / "out_second.csv"
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

    second_exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(second_output),
            "--log-level",
            "INFO",
        ]
    )
    assert second_exit_code == 0
    assert second_output.exists()
    assert output_csv.read_bytes() == second_output.read_bytes()
    second_digest = hashlib.sha256(second_output.read_bytes()).hexdigest()
    assert second_digest == digest

    second_meta_path = second_output.with_name(second_output.name + ".meta.yaml")
    assert second_meta_path.exists()
    second_meta = yaml.safe_load(second_meta_path.read_text())
    assert second_meta["stats"]["output_sha256"] == second_digest
    assert second_meta["schema"] == "PubMedDocumentsSchema"
