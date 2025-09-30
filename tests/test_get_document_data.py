from __future__ import annotations

import argparse
import io
import json
import sys
from collections.abc import Iterable, Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io as lib_io
from library import rate_limiter as rl
from library.cli import LoggerConfig, configure_logger
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
)
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


    def fake_write_csv_chunks(

        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: Sequence[str],
        chunk_size: int | None = None,
        **_: Any,
    ) -> Path:
        list(chunks)
        return path

    def fake_write_export_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: Sequence[str],
        chunk_size: int | None = None,
    ) -> Path:
        list(chunks)
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "_write_export_chunks", fake_write_export_chunks)
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

    cfg = Config()
    chunk_size = cfg.document.all.chunk_size
    ids = [f"CHEMBL{i}" for i in range(1, chunk_size + 3)]
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\n" + "\n".join(ids) + "\n")

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
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
    )
    rc = gdd.run_all(cfg, args)
    assert rc == 1
    log_output = buffer.getvalue()
    configure_logger(LoggerConfig(stream=sys.stdout))
    records = [json.loads(line) for line in log_output.splitlines() if line.strip()]
    assert any(
        record.get("event") == "chembl_documents_fetch_failed"
        and record.get("ids") == ids[:chunk_size]
        for record in records
    )


def test_run_all_passes_generator_to_get_documents(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``run_all`` should pass a generator into ``get_documents``."""

    cfg = Config()
    ids = [f"CHEMBL{i}" for i in range(1, 4)]
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\n" + "\n".join(ids) + "\n")

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    captured: dict[str, object] = {}

    def fake_get_documents(
        ids_iter: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        assert not isinstance(ids_iter, list)
        assert isinstance(ids_iter, Iterator)
        values = list(ids_iter)
        captured["values"] = values
        return pd.DataFrame({"document_chembl_id": values})

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)
    monkeypatch.setattr(gdd.dp, "postprocess_documents", lambda df: df)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_finalise_export(*args: Any, **kwargs: Any) -> int:
        return 0

    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)
    monkeypatch.setattr(
        gdd,
        "fetch_pubmed_records",
        lambda *args, **kwargs: pytest.fail("unexpected call"),
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
    )

    rc = gdd.run_all(cfg, args)
    assert rc == 0
    assert captured["values"] == ids


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
        lambda df, output, cfg, input_csv, key_columns, **kwargs: 0,
    )

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
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


def test_build_fallback_doi_map() -> None:
    """Fallback DOI helper should normalise values and ignore blanks."""

    df = pd.DataFrame(
        {
            "PMID": ["1", "2", "3", None],
            "DOI": ["10.1000/XYZ", "", None, "10.2000/abc"],
        }
    )

    result = gdd._build_fallback_doi_map(df, pmid_column="PMID", doi_column="DOI")

    assert result == {"1": "10.1000/xyz"}


def test_run_pubmed_uses_fallback_csv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PubMed pipeline should forward fallback DOI mappings when provided."""

    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n1\n2\n")
    fallback_csv = tmp_path / "fallback.csv"
    fallback_csv.write_text("PMID,DOI\n1,10.1000/FOO\n2,\n")

    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)
    monkeypatch.setattr(
        gdd,
        "_finalise_export",
        lambda df, output, cfg, input_csv, key_columns, **kwargs: 0,
    )

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
        captured["fallback_doi_map"] = fallback_doi_map
        return pd.DataFrame({"PMID": list(pmids)})

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=fallback_csv,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
    )

    exit_code = gdd.run_pubmed(cfg, args)

    assert exit_code == 0
    assert captured["fallback_doi_map"] == {"1": "10.1000/foo"}


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


    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **_: Any,
    ) -> Path:
        frames = list(chunks)
        captured["col_order"] = list(col_order or [])
        captured["columns"] = list(frames[0].columns) if frames else []
        captured["key_cols"] = list(key_cols or [])
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
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

    factory_calls: dict[str, Any] = {}

    class DummySession:
        def __enter__(self) -> DummySession:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    def fake_session_factory(api_cfg: Any, retry_cfg: Any) -> DummySession:
        factory_calls["api"] = api_cfg
        factory_calls["retry"] = retry_cfg
        return DummySession()

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_factory)

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        client: Any | None = None,
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
    assert factory_calls["api"] is config.api
    assert factory_calls["retry"] is config.retry
    assert "PubMed.PMID" in df.columns
    assert "publication_class" in df.columns
    assert df.loc[0, "PubMed.PMID"] == "1"


def test_fetch_pubmed_records_acquires_documents_limiter(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """PubMed and Semantic Scholar use the shared documents limiter."""

    class TrackingLimiter:
        def __init__(self) -> None:
            self.acquisitions = 0
            self.history: list[int] = []

        def acquire(self) -> None:
            self.acquisitions += 1
            self.history.append(self.acquisitions)

    tracking_limiter = TrackingLimiter()
    batches_seen: list[list[str]] = []
    service_events: list[tuple[str, int]] = []

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd.requests, "Session", lambda: DummySession())

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        batches_seen.append(list(batch))
        service_events.append(("pubmed", tracking_limiter.acquisitions))
        return [{"PubMed.PMID": pmid} for pmid in batch]

    def fake_semantic_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        service_events.append(("semantic", tracking_limiter.acquisitions))
        return [{"scholar.PMID": pmid} for pmid in pmids]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", lambda *_, **__: {})

    def fake_get_limiter(name: str, *_, **__) -> Any:
        if name == "documents_global":
            return tracking_limiter
        return DummyLimiter()

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

    df = gdd.fetch_pubmed_records(
        ["1", "2"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert batches_seen == [["1"], ["2"]]
    assert tracking_limiter.acquisitions == 4
    assert tracking_limiter.history == [1, 2, 3, 4]
    assert service_events == [
        ("pubmed", 1),
        ("semantic", 2),
        ("pubmed", 3),
        ("semantic", 4),
    ]
    assert list(df["PubMed.PMID"]) == ["1", "2"]


def test_documents_limiter_enforces_shared_pace(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Successive service calls respect the shared documents pace."""

    current_time = {"value": 0.0}
    sleep_calls: list[float] = []

    def fake_monotonic() -> float:
        return current_time["value"]

    def fake_sleep(delay: float) -> None:
        sleep_calls.append(delay)
        current_time["value"] += delay

    monkeypatch.setattr(rl.time, "monotonic", fake_monotonic)
    monkeypatch.setattr(rl, "sleep", fake_sleep)

    documents_limiter = rl.RateLimiter(rps=1, burst=1)

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [{"PubMed.PMID": pmid} for pmid in batch]

    def fake_semantic_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [{"scholar.PMID": pmid} for pmid in pmids]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", lambda *_, **__: {})

    def fake_get_limiter(name: str, *_, **__) -> Any:
        if name == "documents_global":
            return documents_limiter
        return DummyLimiter()

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

    df = gdd.fetch_pubmed_records(
        ["1", "2"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert list(df["PubMed.PMID"]) == ["1", "2"]
    assert len(sleep_calls) >= 2
    assert sleep_calls[0] == pytest.approx(1.0, rel=1e-6)
    assert sleep_calls[1] == pytest.approx(1.0, rel=1e-6)

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

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

    seen_cfg: dict[str, PubMedCfg | None] = {"value": None}

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: PubMedCfg | None = None,
        *,
        client: Any | None = None,
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

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        client: Any | None = None,
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


def test_fetch_pubmed_records_falls_back_to_single_semantic_call(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Semantic Scholar batch failures should be retried via the single endpoint."""

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        client: Any | None = None,
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
        return [
            {
                "scholar.PMID": "123",
                "scholar.Error": "Bad Request",
                "scholar.DOI": "",
            }
        ]

    single_calls: list[str] = []

    def fake_semantic_single(
        session: Any, pmid: str, sleep: float, cfg: Any | None = None
    ) -> dict[str, str]:
        single_calls.append(pmid)
        return {
            "scholar.PMID": pmid,
            "scholar.DOI": "10.1000/fallback",
            "scholar.Error": "",
        }

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", fake_semantic_single)
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", lambda *_, **__: {})
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    df = gdd.fetch_pubmed_records(
        ["123"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert single_calls == ["123"]
    assert df.loc[0, "scholar.DOI"] == "10.1000/fallback"


def test_finalise_export_falls_back_to_default_key(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CSV export should default to known key columns when none are provided."""

    cfg = Config()
    df = pd.DataFrame({"document_chembl_id": ["CHEMBL1"], "PubMed.PMID": ["123"]})
    output = tmp_path / "documents.csv"
    captured: dict[str, Any] = {}

    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: list[str] | None = None,
        **_: Any,
    ) -> Path:
        list(chunks)
        captured["key_cols"] = list(key_cols) if key_cols is not None else None
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "hash")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "build_quality_report", lambda df: {})
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gdd.DocumentsSchema, "validate", lambda frame, lazy=True: frame)

    exit_code = gdd._finalise_export(
        df,
        output,
        cfg,
        input_csv=tmp_path / "input.csv",
        key_columns=None,
    )

    assert exit_code == 0
    assert captured["key_cols"] == ["ChEMBL.document_chembl_id"]


def test_finalise_export_skip_quality_reports(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"document_chembl_id": ["CHEMBL1"]})
    output = tmp_path / "documents.csv"
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n")

    class _DummyColumn:
        required = True

    class DummySchema:
        columns = {"document_chembl_id": _DummyColumn()}

        @staticmethod
        def validate(frame: pd.DataFrame, lazy: bool = True) -> pd.DataFrame:
            return frame

    monkeypatch.setattr(gdd, "add_pipeline_metadata", lambda frame: frame)
    monkeypatch.setattr(gdd, "build_dataframe", lambda frame, **__: frame)
    monkeypatch.setattr(gdd, "DocumentsSchema", DummySchema)
    monkeypatch.setattr(gdd, "_EXPORT_COLUMNS", ["document_chembl_id"])
    monkeypatch.setattr(gdd, "_EXPORT_COLUMN_RENAMES", {})
    monkeypatch.setattr(gdd, "_EXPORT_SORT_FALLBACK", ["document_chembl_id"])
    monkeypatch.setattr(gdd, "_EXPORT_STREAM_CHUNK_SIZE", 1)
    monkeypatch.setattr(gdd, "_iter_export_chunks", lambda frame, chunk_size: (frame,))

    def fake_write_chunks(
        chunks: Iterable[pd.DataFrame], path: Path, *, cfg: Any, **__: Any
    ) -> Path:
        frame = next(iter(chunks))
        path = Path(path)
        frame.to_csv(path, index=False)
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda _: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)

    rc = gdd._finalise_export(
        df,
        output,
        cfg,
        input_csv=input_csv,
        key_columns=["document_chembl_id"],
        quality_enabled=False,
    )

    assert rc == 0
    quality_json = output.with_suffix(".quality.json")
    base = output.with_suffix("")
    quality_tables = [
        base.with_name(f"{base.name}_quality_report_table.csv"),
        base.with_name(f"{base.name}_data_correlation_report_table.csv"),
    ]
    for artefact in [quality_json, *quality_tables]:
        assert not artefact.exists()


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
        session, batch, sleep, cfg=None, *, client=None
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
    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

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
