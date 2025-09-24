from __future__ import annotations

import argparse
import io
import sys
import time
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from contextlib import contextmanager

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io as lib_io
from library.cli import LoggerConfig, configure_logger
from library.config import Config
from schemas import DocumentsSchema
from scripts import get_document_data as gdd


def test_cli_uses_custom_column(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure CLI respects the ``--column`` option for ChEMBL documents."""
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n")

    captured: dict[str, str] = {}
    orig_read_ids = lib_io.read_ids

    def fake_read_ids(
        path: str | Path,
        *,
        column: str,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Iterable[str]:
        captured["column"] = column
        return orig_read_ids(path, column=column, cfg=cfg, sep=sep, encoding=encoding)

    monkeypatch.setattr(lib_io, "read_ids", fake_read_ids)

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    def fake_get_documents(
        ids: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        return pd.DataFrame({"document_chembl_id": data, "title": ["t"] * len(data)})

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_write_csv(
        df: pd.DataFrame,
        path: Path,
        *,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        chunksize: int | None = None,
    ) -> Path:
        return path

    monkeypatch.setattr(lib_io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gdd, "ensure_dirs", lambda cfg: None)

    rc = gdd.main(
        [
            "chembl",
            "--column",
            "document_chembl_id",
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
        ]
    )
    assert rc == 0
    assert captured["column"] == "document_chembl_id"


def test_run_all_logs_failing_ids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Log message should include IDs when document retrieval fails."""
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\nCHEMBL2\n")

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    def fake_get_documents(
        ids: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        raise ValueError("boom")

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)

    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="ERROR", stream=buffer))
    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
    )
    rc = gdd.run_all(cfg, args)
    assert rc == 1
    log_output = buffer.getvalue()
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert "CHEMBL1" in log_output and "CHEMBL2" in log_output


def test_write_csv_column_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Columns follow schema order with extras alphabetically appended."""
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n")

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    df = pd.DataFrame(
        {
            "title": ["t"],
            "B": ["2"],
            "document_chembl_id": ["CHEMBL1"],
            "A": ["1"],
            "doi": ["10.0"],
        }
    )

    monkeypatch.setattr(cl, "get_documents", lambda *_, **__: df)
    monkeypatch.setattr(gdd, "normalize_documents", lambda frame: frame)

    captured: dict[str, list[str]] = {}

    def fake_write_csv(
        frame: pd.DataFrame,
        path: Path,
        *,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        chunksize: int | None = None,
    ) -> Path:
        captured["col_order"] = list(col_order or [])
        return path

    monkeypatch.setattr(lib_io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
    )
    rc = gdd.run_chembl(cfg, args)
    assert rc == 0
    schema_cols = list(DocumentsSchema.columns)
    expected = [c for c in schema_cols if c in df.columns] + sorted(
        c for c in df.columns if c not in schema_cols
    )
    assert captured["col_order"] == expected


def test_fetch_pubmed_records_order_and_limiters(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    pmids = ["1", "2"]

    acquisitions: list[str] = []
    limiters: dict[str, FakeLimiter] = {}

    class FakeLimiter:
        def __init__(self, name: str) -> None:
            self.name = name

        def acquire(self) -> None:
            acquisitions.append(self.name)

    def fake_get_limiter(name: str, rps: float, burst: int | None = None) -> FakeLimiter:
        limiter = limiters.get(name)
        if limiter is None:
            limiter = FakeLimiter(name)
            limiters[name] = limiter
        return limiter

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

    def fake_fetch_pubmed_batch(
        session: object,
        batch: list[str],
        sleep: float,
        cfg: object | None = None,
    ) -> list[dict[str, str]]:
        pmid = batch[0]
        if pmid == "1":
            time.sleep(0.02)
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"doi-{pmid}",
                "PubMed.Error": "",
            }
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)

    def fake_semantic_batch(
        session: object,
        batch_pmids: list[str],
        sleep: float,
        cfg: object | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"doi-{pmid}",
                "scholar.Error": "",
            }
            for pmid in batch_pmids
        ]

    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar_batch",
        fake_semantic_batch,
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda session, pmid, cfg, limiter: (
            limiter.acquire(),
            {"OpenAlex.Id": f"oa-{pmid}"},
        )[1],
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda session, doi, cfg, limiter: (
            limiter.acquire(),
            {"crossref.Type": "journal"},
        )[1],
    )

    @contextmanager
    def fake_session_with_retry(api_cfg: object, retry_cfg: object) -> object:
        yield object()

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_with_retry)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        openalex_cfg=cfg.openalex,
        crossref_cfg=cfg.crossref,
        api_cfg=cfg.api,
        retry_cfg=cfg.retry,
        pubmed_cfg=cfg.pubmed,
        semantic_cfg=cfg.semantic_scholar,
        max_workers=2,
        batch_size=1,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert acquisitions.count("pubmed") == len(pmids)
    assert acquisitions.count("semantic_scholar") == len(pmids)
    assert acquisitions.count("openalex") == len(pmids)
    assert acquisitions.count("crossref") == len(pmids)
    assert set(limiters) == {"pubmed", "semantic_scholar", "openalex", "crossref"}
