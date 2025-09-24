from __future__ import annotations

import argparse
import io
import json
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io as lib_io
from library.cli import LoggerConfig, configure_logger
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
)
from library.document_pipeline import DOCUMENT_SCHEMA_COLUMNS
from schemas import DocumentsSchema
from scripts import get_document_data as gdd


class DummyLimiter:
    def acquire(self) -> None:
        return None


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
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "build_quality_report", lambda df: {})
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
    records = [json.loads(line) for line in log_output.splitlines() if line.strip()]
    assert any(
        record.get("event") == "chembl_documents_fetch_failed"
        and record.get("ids") == ["CHEMBL1", "CHEMBL2"]
        for record in records
    )


def test_pubmed_cli_rejects_non_positive_batch_size(tmp_path: Path) -> None:
    """PubMed CLI should reject zero or negative batch sizes."""

    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n1\n")

    with pytest.raises(SystemExit) as exc_info:
        gdd.main(
            [
                "pubmed",
                "--input",
                str(input_csv),
                "--output",
                str(tmp_path / "out.csv"),
                "--batch-size",
                "0",
            ]
        )
    assert exc_info.value.code == 2


def test_chembl_cli_rejects_non_positive_chunk_size(tmp_path: Path) -> None:
    """ChEMBL CLI should reject zero or negative chunk sizes."""

    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n")

    with pytest.raises(SystemExit) as exc_info:
        gdd.main(
            [
                "chembl",
                "--input",
                str(input_csv),
                "--output",
                str(tmp_path / "out.csv"),
                "--chunk-size",
                "0",
            ]
        )
    assert exc_info.value.code == 2


def test_run_pubmed_uses_keyword_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PubMed sub-command forwards configuration via keyword arguments."""

    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n1\n2\n")

    def fake_read_ids(
        path: str | Path,
        *,
        column: str,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
        na_markers: list[str] | None = None,
    ) -> Iterable[str]:
        assert Path(path) == input_csv
        assert column == "PMID"
        return iter(["1", "2"])

    monkeypatch.setattr(lib_io, "read_ids", fake_read_ids)

    captured: dict[str, Any] = {}

    def fake_fetch_pubmed_records(
        pmids: Iterable[str],
        cfg_param: Config,
        *,
        sleep: float,
        pubmed_cfg: Any | None = None,
        semantic_scholar_cfg: SemanticScholarCfg,
        openalex_cfg: OpenAlexCfg,
        crossref_cfg: CrossRefCfg,
        max_workers: int,
        batch_size: int,
        fallback_doi_map: dict[str, str] | None = None,
    ) -> pd.DataFrame:
        captured["pmids"] = list(pmids)
        captured["cfg"] = cfg_param
        captured["sleep"] = sleep
        captured["pubmed_cfg"] = pubmed_cfg
        captured["semantic_scholar_cfg"] = semantic_scholar_cfg
        captured["openalex_cfg"] = openalex_cfg
        captured["crossref_cfg"] = crossref_cfg
        captured["max_workers"] = max_workers
        captured["batch_size"] = batch_size
        captured["fallback_doi_map"] = fallback_doi_map
        return pd.DataFrame({"PMID": ["1", "2"]})

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)
    monkeypatch.setattr(
        gdd,
        "_finalise_export",
        lambda df, output, cfg, input_csv, key_columns: 0,
    )

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
    )

    exit_code = gdd.run_pubmed(cfg, args)
    assert exit_code == 0
    assert captured["pmids"] == ["1", "2"]
    assert captured["cfg"] is cfg
    assert captured["sleep"] == cfg.document.pubmed.sleep
    assert captured["pubmed_cfg"] is cfg.pubmed
    assert captured["semantic_scholar_cfg"] is cfg.semantic_scholar
    assert captured["openalex_cfg"] is cfg.openalex
    assert captured["crossref_cfg"] is cfg.crossref
    assert captured["max_workers"] == cfg.document.pubmed.workers
    assert captured["batch_size"] == cfg.document.pubmed.batch_size
    assert captured["fallback_doi_map"] is None


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

    captured: dict[str, Any] = {}

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
        captured["columns"] = list(frame.columns)
        captured["key_cols"] = list(key_cols or [])
        return path

    monkeypatch.setattr(lib_io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "build_quality_report", lambda df: {})

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
    )
    rc = gdd.run_chembl(cfg, args)
    assert rc == 0
    assert captured["col_order"] == gdd._EXPORT_COLUMNS
    assert captured["columns"] == gdd._EXPORT_COLUMNS
    assert captured["key_cols"] == ["ChEMBL.document_chembl_id"]


def test_fetch_pubmed_records_handles_generic_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure unexpected errors yield failure records rather than crashing."""

    def fail_fetch_pubmed_batch(
        *args: object, **kwargs: object
    ) -> list[dict[str, str]]:
        raise RuntimeError("boom")

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fail_fetch_pubmed_batch)
    monkeypatch.setattr(gdd, "get_limiter", lambda *args, **kwargs: DummyLimiter())

    df = gdd.fetch_pubmed_records(
        ["123"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert len(df) == 1
    row = df.iloc[0]
    assert row["PubMed.PMID"] == "123"
    assert row["PubMed.Error"] == "boom"
    assert row["scholar.Error"] == "boom"
    assert row["OpenAlex.Error"] == "boom"
    assert row["crossref.Error"] == "boom"
    assert row["publication_class"] == "unknown"


def test_fetch_pubmed_records_accepts_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Legacy two-argument usage should still be supported."""

    config = Config()
    expected_sleep = config.document.pubmed.sleep
    semantic_cfg = config.semantic_scholar

    def fake_session() -> Any:
        class DummySession:
            def __enter__(self) -> DummySession:
                return self

            def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
                return None

        return DummySession()

    monkeypatch.setattr(gdd.requests, "Session", fake_session)

    def fake_pubmed_batch(
        session: Any, batch: list[str], sleep: float, cfg: Any | None = None
    ) -> list[dict[str, str]]:
        assert sleep == expected_sleep
        assert batch == ["1"]
        assert cfg is config.pubmed
        return [
            {
                "PubMed.PMID": "1",
                "PubMed.DOI": "10.1000/xyz",
                "PubMed.PublicationType": "Review",
            }
        ]

    def fake_semantic_scholar_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        assert pmids == ["1"]
        assert cfg is semantic_cfg
        return [
            {
                "scholar.PMID": "1",
                "scholar.PublicationTypes": "Review",
                "scholar.DOI": "10.1000/xyz",
            }
        ]

    def fake_openalex(
        session: Any, pmid: str, cfg_arg: Any, limiter: Any
    ) -> dict[str, str]:
        assert pmid == "1"
        assert cfg_arg is config.openalex
        return {"OpenAlex.PublicationTypes": "journal-article"}

    def fake_crossref(
        session: Any, doi: str, cfg_arg: Any, limiter: Any
    ) -> dict[str, str]:
        assert doi == "10.1000/xyz"
        assert cfg_arg is config.crossref
        return {"crossref.Type": "journal-article"}

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(
        gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_scholar_batch
    )
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)
    monkeypatch.setattr(gdd, "get_limiter", lambda *args, **kwargs: DummyLimiter())

    df = gdd.fetch_pubmed_records(["1"], config)
    assert "PubMed.PMID" in df.columns
    assert "publication_class" in df.columns
    assert df.loc[0, "PubMed.PMID"] == "1"


def test_fetch_pubmed_records_uses_explicit_pubmed_cfg(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Explicit PubMed configuration should be passed to the batch fetcher."""

    sentinel_cfg = PubMedCfg(base="https://example.org/api")

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd.requests, "Session", lambda: DummySession())

    seen_cfg: dict[str, PubMedCfg | None] = {"value": None}

    def fake_pubmed_batch(
        session: Any, batch: list[str], sleep: float, cfg: PubMedCfg | None = None
    ) -> list[dict[str, str]]:
        seen_cfg["value"] = cfg
        return [{"PubMed.PMID": pmid} for pmid in batch]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", lambda *_, **__: [])
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", lambda *_, **__: {})
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    df = gdd.fetch_pubmed_records(
        ["1"],
        sleep=0.0,
        pubmed_cfg=sentinel_cfg,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert seen_cfg["value"] is sentinel_cfg
    assert df.loc[0, "PubMed.PMID"] == "1"


def test_fetch_pubmed_records_uses_fallback_doi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CrossRef queries should fall back to supplied DOI mapping when needed."""

    fallback = {"123": "10.9999/fallback"}

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd.requests, "Session", lambda: DummySession())

    def fake_pubmed_batch(
        session: Any, batch: list[str], sleep: float, cfg: Any | None = None
    ) -> list[dict[str, str]]:
        assert batch == ["123"]
        return [{"PubMed.PMID": "123"}]

    def fake_semantic_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        assert pmids == ["123"]
        return [{"scholar.PMID": "123"}]

    def fake_openalex(
        session: Any, pmid: str, cfg_arg: Any, limiter: Any
    ) -> dict[str, str]:
        assert pmid == "123"
        return {}

    captured: dict[str, str] = {}

    def fake_crossref(
        session: Any, doi: str, cfg_arg: Any, limiter: Any
    ) -> dict[str, str]:
        captured["doi"] = doi
        return {"crossref.DOI": doi}

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    df = gdd.fetch_pubmed_records(
        ["123"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
        fallback_doi_map=fallback,
    )

    assert captured["doi"] == fallback["123"]
    assert df.loc[0, "crossref.DOI"] == fallback["123"]


@pytest.mark.parametrize("context_position", ["suffix", "prefix"])
def test_fetch_pubmed_records_accepts_executor_context(
    monkeypatch: pytest.MonkeyPatch,
    context_position: str,
) -> None:
    """Executor passing an internal context argument should be ignored."""

    pmids = ["1", "2"]

    class DummyLimiter:
        def acquire(self) -> None:  # pragma: no cover - trivial
            return None

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    class DummyFuture:
        def __init__(self, value: list[dict[str, str]]) -> None:
            self._value = value

        def result(self) -> list[dict[str, str]]:  # pragma: no cover - trivial
            return self._value

        def __hash__(self) -> int:  # pragma: no cover - deterministic identity
            return id(self)

    class DummyExecutor:
        def __init__(self, *args: object, **kwargs: object) -> None:
            self._submitted: list[DummyFuture] = []

        def __enter__(self) -> DummyExecutor:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

        def submit(self, fn, batch):  # type: ignore[no-untyped-def]
            context = object()
            if context_position == "prefix":
                value = fn(context, batch)
            else:
                value = fn(batch, context)
            future = DummyFuture(value)

            self._submitted.append(future)
            return future

    def fake_build_dataframe(records, **_):  # type: ignore[no-untyped-def]
        return pd.DataFrame.from_records(records)

    def fake_merge_metadata(pubmed, semsch, openalex, crossref):  # type: ignore[no-untyped-def]
        combined: dict[str, str] = {}
        combined.update(pubmed)
        combined.update(semsch)
        combined.update(openalex)
        combined.update(crossref)
        return combined

    def fake_pubmed_batch(  # type: ignore[no-untyped-def]
        session, batch, sleep, cfg=None
    ):
        return [{"PubMed.PMID": pmid, "PubMed.DOI": pmid} for pmid in batch]

    def fake_sem_batch(session, batch, sleep, cfg=None):  # type: ignore[no-untyped-def]
        return [{"scholar.PMID": pmid} for pmid in batch]

    def fake_openalex(session, pmid, cfg, limiter):  # type: ignore[no-untyped-def]
        return {"OpenAlex.ID": pmid}

    def fake_crossref(session, doi, cfg, limiter):  # type: ignore[no-untyped-def]
        return {"crossref.DOI": doi}

    monkeypatch.setattr(gdd, "ThreadPoolExecutor", DummyExecutor)
    monkeypatch.setattr(gdd, "as_completed", lambda futures: iter(futures))
    monkeypatch.setattr(gdd, "build_dataframe", fake_build_dataframe)
    monkeypatch.setattr(gdd, "merge_metadata", fake_merge_metadata)
    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_sem_batch)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())
    monkeypatch.setattr(gdd.requests, "Session", DummySession)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=2,
    )

    assert sorted(df["PubMed.PMID"].tolist()) == pmids
