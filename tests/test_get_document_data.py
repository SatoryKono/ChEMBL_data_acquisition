from __future__ import annotations

import argparse
import gc
import io
import shutil
import json
import sys
import threading
import time
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from concurrent.futures import Future
from pathlib import Path
from typing import Any, Callable

from contextlib import contextmanager

from itertools import islice

import pandas as pd
import pytest
import requests

responses = pytest.importorskip("responses")

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
    def __init__(self, burst: int = 1) -> None:
        self.burst = burst

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
        frames = list(chunks)
        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame()
        combined.to_csv(path, index=False)
        return path

    def fake_write_export_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: Sequence[str],
        chunk_size: int | None = None,
    ) -> Path:
        frames = list(chunks)
        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame()
        combined.to_csv(path, index=False)
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "_write_export_chunks", fake_write_export_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(
        gdd, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "build_quality_report", lambda *_, **__: {})
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


def test_run_pubmed_large_limit_streams(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Large limits should not force PubMed runs to materialise identifiers."""

    cfg = Config()
    limit = 1_000_000
    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n1\n")

    def fake_read_ids(
        *_: Any,
        **__: Any,
    ) -> Iterable[str]:
        return (str(i) for i in range(limit))

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
        fallback_doi_map: Mapping[str, str] | None = None,
        return_generator: bool = False,
    ) -> Iterable[pd.DataFrame]:
        assert not isinstance(pmids, list)
        assert isinstance(pmids, Iterator)
        subset = list(islice(pmids, 3))
        captured["pmids"] = subset

        def _generator() -> Iterator[pd.DataFrame]:
            yield pd.DataFrame({"PMID": subset})

        return _generator()

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_finalise_export(
        frames: Iterable[pd.DataFrame],
        output: Path,
        cfg_param: Config,
        *,
        input_csv: Path,
        key_columns: Sequence[str] | None,
        **kwargs: Any,
    ) -> int:
        collected = list(frames)
        assert len(collected) == 1
        assert collected[0]["PMID"].tolist() == captured["pmids"]
        return 0

    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)

    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="INFO", stream=buffer))

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
        limit=limit,
    )

    exit_code = gdd.run_pubmed(cfg, args)
    assert exit_code == 0
    assert captured["pmids"] == ["0", "1", "2"]

    records = [
        json.loads(line)
        for line in buffer.getvalue().splitlines()
        if line.strip()
    ]
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert any(
        record.get("event") == "process_limit" and record.get("limit") == 3
        for record in records
    )


def test_run_all_large_limit_streams(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Large limits should leave ``run_all`` streaming identifiers lazily."""

    cfg = Config()
    limit = 1_000_000
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL0\n")

    def fake_read_ids(
        *_: Any,
        **__: Any,
    ) -> Iterable[str]:
        return (f"CHEMBL{i}" for i in range(limit))

    monkeypatch.setattr(lib_io, "read_ids", fake_read_ids)

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - no cleanup
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    captured: dict[str, Any] = {}

    def fake_get_documents(
        ids_iter: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        assert not isinstance(ids_iter, list)
        assert isinstance(ids_iter, Iterator)
        values = list(islice(ids_iter, 5))
        captured["values"] = values
        return pd.DataFrame(
            {
                "document_chembl_id": values,
                "pubmed_id": list(range(1, 6)),
                "doi": [f"10.1000/{i}" for i in range(1, 6)],
            }
        )

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)
    monkeypatch.setattr(gdd, "merge_with_chembl", lambda doc_df, _: doc_df)
    monkeypatch.setattr(gdd.dp, "postprocess_documents", lambda df: df)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)
    monkeypatch.setattr(gdd, "fetch_pubmed_records", lambda *args, **kwargs: iter(()))

    def fake_finalise_export(*_: Any, **__: Any) -> int:
        return 0

    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)

    buffer = io.StringIO()
    configure_logger(LoggerConfig(level="INFO", stream=buffer))

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
        limit=limit,
        chunk_size=5,
    )

    exit_code = gdd.run_all(cfg, args)
    assert exit_code == 0
    assert captured["values"] == [f"CHEMBL{i}" for i in range(5)]

    records = [
        json.loads(line)
        for line in buffer.getvalue().splitlines()
        if line.strip()
    ]
    configure_logger(LoggerConfig(stream=sys.stdout))
    assert any(
        record.get("event") == "process_limit" and record.get("limit") == 5
        for record in records
    )


def test_run_all_streams_pubmed_batches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PubMed batches are written to disk and merged chunk by chunk."""

    cfg = Config()
    input_csv = tmp_path / "docs.csv"
    ids = [f"CHEMBL{i}" for i in range(40)]
    input_csv.write_text("document_chembl_id\n" + "\n".join(ids))

    monkeypatch.setattr(lib_io, "read_ids", lambda *_args, **_kwargs: iter(ids))

    class DummyClient:
        def __enter__(self) -> "DummyClient":  # pragma: no cover - simple context
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - simple context
            return None

    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: DummyClient())

    def fake_get_documents(
        ids_iter: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        values = list(ids_iter)
        return pd.DataFrame(
            {
                "document_chembl_id": values,
                "pubmed_id": list(range(len(values))),
                "doi": [f"10.1000/{idx}" for idx in range(len(values))],
            }
        )

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)

    batch_size = 8

    def fake_fetch_pubmed_records(
        pmids: Sequence[str], *args: object, **kwargs: object
    ) -> Iterator[pd.DataFrame]:
        def _generator() -> Iterator[pd.DataFrame]:
            for start in range(0, len(pmids), batch_size):
                chunk_pmids = pmids[start : start + batch_size]
                yield pd.DataFrame(
                    {
                        "PubMed.PMID": chunk_pmids,
                        "PubMed.DOI": [f"10.2000/{pmid}" for pmid in chunk_pmids],
                    }
                )

        return _generator()

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)

    cleanup_flag = {"cleaned": False}

    def tracking_tempdir(*_args: object, **_kwargs: object):
        base = tmp_path / "pubmed_tmp"
        if base.exists():
            shutil.rmtree(base)
        base.mkdir()

        class _Context:
            def __enter__(self) -> str:  # pragma: no cover - simple context
                return str(base)

            def __exit__(
                self, exc_type: object, exc: object, tb: object | None
            ) -> None:  # pragma: no cover - simple context
                shutil.rmtree(base)
                cleanup_flag["cleaned"] = True

        return _Context()

    monkeypatch.setattr(gdd.tempfile, "TemporaryDirectory", tracking_tempdir)

    metadata_arguments: list[object] = []
    real_merge = gdd.merge_with_chembl

    def recording_merge(doc_df: pd.DataFrame, metadata: object) -> pd.DataFrame:
        metadata_arguments.append(metadata)
        return real_merge(doc_df, metadata)

    monkeypatch.setattr(gdd, "merge_with_chembl", recording_merge)

    captured: dict[str, pd.DataFrame] = {}

    def fake_postprocess(df: pd.DataFrame) -> pd.DataFrame:
        captured["merged"] = df.copy()
        return df

    monkeypatch.setattr(gdd.dp, "postprocess_documents", fake_postprocess)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_finalise_export(
        df: pd.DataFrame | Iterable[pd.DataFrame], *args: object, **kwargs: object
    ) -> int:
        assert isinstance(df, pd.DataFrame)
        captured["export_rows"] = len(df)
        return 0

    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)

    monkeypatch.setattr(
        gdd.pd,
        "concat",
        lambda *_args, **_kwargs: pytest.fail("pd.concat not expected"),
    )

    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
        fallback_doi_csv=None,
        fallback_doi_pmid_column="PMID",
        fallback_doi_value_column="DOI",
        limit=None,
        offset=0,
        chunk_size=batch_size,
        batch_size=batch_size,
        column="document_chembl_id",
    )

    exit_code = gdd.run_all(cfg, args)
    assert exit_code == 0
    assert cleanup_flag["cleaned"] is True
    assert metadata_arguments, "merge_with_chembl should receive streamed metadata"
    metadata_obj = metadata_arguments[0]
    assert not isinstance(metadata_obj, pd.DataFrame)
    assert hasattr(metadata_obj, "__iter__")
    merged = captured["merged"]
    assert len(merged) == len(ids)
    assert set(merged["PubMed.DOI"]) == {
        f"10.2000/{idx}" for idx in map(str, range(len(ids)))
    }
    assert captured["export_rows"] == len(ids)


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
        return_generator: bool = False,
    ) -> Iterable[pd.DataFrame]:
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
        assert return_generator is True

        def _generator() -> Iterator[pd.DataFrame]:
            yield pd.DataFrame({"PMID": ["1", "2"]})

        return _generator()

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_finalise_export(
        frames: Iterable[pd.DataFrame],
        output: Path,
        cfg: Config,
        *,
        input_csv: Path,
        key_columns: Sequence[str] | None,
        **kwargs: Any,
    ) -> int:
        assert not isinstance(frames, pd.DataFrame)
        collected = list(frames)
        assert len(collected) == 1
        assert collected[0]["PMID"].tolist() == ["1", "2"]
        return 0

    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)

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
        return_generator: bool = False,
    ) -> Iterable[pd.DataFrame]:
        captured["fallback_doi_map"] = fallback_doi_map
        assert return_generator is True

        def _generator() -> Iterator[pd.DataFrame]:
            yield pd.DataFrame({"PMID": list(pmids)})

        return _generator()

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


def test_run_pubmed_streams_batches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PubMed CLI should stream chunks to the exporter."""

    input_csv = tmp_path / "pmids.csv"
    input_csv.write_text("PMID\n1\n2\n3\n4\n5\n6\n")

    monkeypatch.setattr(
        lib_io,
        "read_ids",
        lambda *_, **__: iter(["1", "2", "3", "4", "5", "6"]),
    )

    chunk_lengths: list[int] = []

    def fake_fetch_pubmed_records(
        *_,
        return_generator: bool = False,
        **__: Any,
    ) -> Iterable[pd.DataFrame]:
        assert return_generator is True

        def _generator() -> Iterator[pd.DataFrame]:
            for start in range(0, 6, 2):
                yield pd.DataFrame({"PMID": [str(start + 1), str(start + 2)]})

        return _generator()

    def fake_normalize(df: pd.DataFrame) -> pd.DataFrame:
        chunk_lengths.append(len(df))
        return df

    def fake_finalise_export(
        frames: Iterable[pd.DataFrame],
        output: Path,
        cfg: Config,
        *,
        input_csv: Path,
        key_columns: Sequence[str] | None,
        **kwargs: Any,
    ) -> int:
        lengths = [len(frame) for frame in frames]
        assert lengths == [2, 2, 2]
        return 0

    monkeypatch.setattr(gdd, "fetch_pubmed_records", fake_fetch_pubmed_records)
    monkeypatch.setattr(gdd, "normalize_documents", fake_normalize)
    monkeypatch.setattr(gdd, "_finalise_export", fake_finalise_export)

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
    assert chunk_lengths == [2, 2, 2]


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
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        combined.to_csv(path, index=False)
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(
        gdd, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "build_quality_report", lambda *_, **__: {})

    cfg = Config()
    args = argparse.Namespace(
        input_csv=input_csv,
        output_csv=tmp_path / "out.csv",
    )
    rc = gdd.run_chembl(cfg, args)
    assert rc == 0
    assert captured["columns"] == gdd._EXPORT_COLUMNS
    assert captured["key_cols"] == ["ChEMBL.document_chembl_id"]


def test_prepare_export_frame_merges_duplicate_columns() -> None:
    """Export preparation should not emit duplicate ChEMBL columns."""

    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "ChEMBL.title": ["", "existing"],
            "title": ["fallback", "ignored"],
            "ChEMBL.abstract": ["pref", ""],
            "abstract": ["", "fallback abstract"],
        }
    )

    result = gdd._prepare_export_frame(df)

    duplicates = [name for name, count in Counter(result.columns).items() if count > 1]
    assert duplicates == []
    assert result.loc[0, "ChEMBL.title"] == "fallback"
    assert result.loc[1, "ChEMBL.title"] == "existing"
    assert result.loc[0, "ChEMBL.abstract"] == "pref"
    assert result.loc[1, "ChEMBL.abstract"] == "fallback abstract"


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


def test_fetch_pubmed_records_returns_generator_in_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Generator mode should yield batches in submission order."""

    pmids = ["1", "2", "3", "4"]

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
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"PubMed.PMID": pmid, "PubMed.DOI": f"10.1000/{pmid}"} for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        ids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"scholar.PMID": pmid, "scholar.DOI": f"10.1000/{pmid}"} for pmid in ids
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_, **__: {"scholar.DOI": "10.1000/fallback"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda session, pmid, cfg, limiter: {"OpenAlex.Id": f"OA{pmid}"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda session, doi, cfg, limiter: {"crossref.DOI": doi},
    )
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    generator = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=2,
        return_generator=True,
    )

    frames = list(generator)
    assert [frame["PubMed.PMID"].tolist() for frame in frames] == [
        ["1", "2"],
        ["3", "4"],
    ]

    combined = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=2,
    )

    assert combined["PubMed.PMID"].tolist() == pmids


def test_fetch_pubmed_records_closes_sessions_when_generator_gc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Generator GC should release thread resources and close sessions."""

    pmids = ["1", "2", "3"]

    class TrackingSession:
        def __init__(self) -> None:
            self.closed = False

        def close(self) -> None:
            self.closed = True

    class TrackingContext:
        def __init__(self, session: TrackingSession) -> None:
            self._session = session

        def __enter__(self) -> TrackingSession:
            sessions.append(self._session)
            return self._session

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - deterministic
            self._session.close()

    sessions: list[TrackingSession] = []

    monkeypatch.setattr(
        gdd,
        "session_with_retry",
        lambda *_, **__: TrackingContext(TrackingSession()),
    )

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"PubMed.PMID": pmid, "PubMed.DOI": f"10.1000/{pmid}"} for pmid in batch
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar_batch",
        lambda *_, **__: [],
    )
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    generator = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
        return_generator=True,
    )

    next(generator)
    assert any(not session.closed for session in sessions)

    del generator
    gc.collect()

    deadline = time.time() + 1.0
    while any(not session.closed for session in sessions) and time.time() < deadline:
        gc.collect()
        time.sleep(0.01)

    assert sessions and all(session.closed for session in sessions)


def test_fetch_pubmed_records_drains_pending_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """More batches than ``max_in_flight`` should still all be emitted."""

    pmids = [str(i) for i in range(6)]

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
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"PubMed.PMID": pmid, "PubMed.DOI": f"10.1000/{pmid}"} for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        ids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"scholar.PMID": pmid, "scholar.DOI": f"10.1000/{pmid}"} for pmid in ids
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_, **__: {"scholar.DOI": "10.1000/fallback"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda session, pmid, cfg, limiter: {"OpenAlex.Id": f"OA{pmid}"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda session, doi, cfg, limiter: {"crossref.DOI": doi},
    )
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    generator = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=2,
        return_generator=True,
    )

    frames = list(generator)
    assert [frame["PubMed.PMID"].tolist() for frame in frames] == [
        ["0", "1"],
        ["2", "3"],
        ["4", "5"],
    ]

    combined = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=2,
    )

    assert combined["PubMed.PMID"].tolist() == pmids


def test_fetch_pubmed_records_limits_heap_with_delayed_first_batch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Completed results stay bounded when the first batch is slow."""

    pmids = [str(i) for i in range(40)]
    max_workers = 2

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    first_batch_started = threading.Event()
    other_batch_ready = threading.Event()
    release_first_batch = threading.Event()

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        if batch and batch[0] == "0":
            first_batch_started.set()
            release_first_batch.wait()
        return [
            {"PubMed.PMID": pmid, "PubMed.DOI": f"10.1000/{pmid}"} for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        ids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"scholar.PMID": pmid, "scholar.DOI": f"10.1000/{pmid}"} for pmid in ids
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_, **__: {"scholar.DOI": "10.1000/fallback"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda session, pmid, cfg, limiter: {"OpenAlex.Id": f"OA{pmid}"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda session, doi, cfg, limiter: {"crossref.DOI": doi},
    )

    heap_lock = threading.Lock()
    max_heap_size = 0
    heap_sizes_after_release: list[int] = []

    original_heappush = gdd.heapq.heappush
    original_heappop = gdd.heapq.heappop

    def tracking_heappush(
        heap: list[tuple[int, list[dict[str, str]]]],
        item: tuple[int, list[dict[str, str]]],
    ) -> None:
        nonlocal max_heap_size
        original_heappush(heap, item)
        with heap_lock:
            size = len(heap)
            if item[0] != 0 and not release_first_batch.is_set():
                other_batch_ready.set()
            if size > max_heap_size:
                max_heap_size = size
            if release_first_batch.is_set():
                heap_sizes_after_release.append(size)

    def tracking_heappop(
        heap: list[tuple[int, list[dict[str, str]]]],
    ) -> tuple[int, list[dict[str, str]]]:
        item = original_heappop(heap)
        with heap_lock:
            if release_first_batch.is_set():
                heap_sizes_after_release.append(len(heap))
        return item

    monkeypatch.setattr(gdd.heapq, "heappush", tracking_heappush)
    monkeypatch.setattr(gdd.heapq, "heappop", tracking_heappop)

    generator = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=max_workers,
        batch_size=1,
        return_generator=True,
    )

    frames: list[pd.DataFrame] = []

    def consume() -> None:
        for frame in generator:
            frames.append(frame)

    consumer = threading.Thread(target=consume)
    consumer.start()

    try:
        assert first_batch_started.wait(timeout=5)
        assert other_batch_ready.wait(timeout=5)
        time.sleep(0.2)
        release_first_batch.set()
        consumer.join(timeout=5)
    finally:
        release_first_batch.set()
        consumer.join(timeout=5)

    assert not consumer.is_alive()
    assert frames
    combined = pd.concat(frames, ignore_index=True)
    assert combined["PubMed.PMID"].tolist() == pmids

    assert max_heap_size > max_workers
    assert max_heap_size <= max_workers * 4
    assert heap_sizes_after_release and heap_sizes_after_release[-1] == 0


def test_fetch_pubmed_records_streams_large_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Large input should be streamed in bounded chunks."""

    pmids = [str(i) for i in range(30)]

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    def fake_pubmed_batch(
        session: Any,
        batch: list[str],
        sleep: float,
        cfg: Any | None = None,
        *,
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"PubMed.PMID": pmid, "PubMed.DOI": f"10.1000/{pmid}"} for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        ids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {"scholar.PMID": pmid, "scholar.DOI": f"10.1000/{pmid}"} for pmid in ids
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_, **__: {"scholar.DOI": "10.1000/fallback"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda session, pmid, cfg, limiter: {"OpenAlex.Id": f"OA{pmid}"},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda session, doi, cfg, limiter: {"crossref.DOI": doi},
    )

    generator = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=5,
        return_generator=True,
    )

    batch_sizes = [len(frame) for frame in generator]
    assert batch_sizes == [5, 5, 5, 5, 5, 5]


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
        retry_cfg: Any | None = None,
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
        retry_cfg: Any | None = None,
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


def test_openalex_and_crossref_jobs_acquire_service_limiters(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OpenAlex and CrossRef jobs use their service-specific limiters."""

    class TrackingLimiter:
        def __init__(self, name: str) -> None:
            self.name = name
            self.acquisitions = 0
            self._lock = threading.Lock()

        def acquire(self) -> None:
            with self._lock:
                self.acquisitions += 1

    limiters = {
        "documents_global": TrackingLimiter("global"),
        "documents_openalex": TrackingLimiter("openalex"),
        "documents_crossref": TrackingLimiter("crossref"),
    }

    def fake_get_limiter(name: str, *_, **__) -> Any:
        return limiters.get(name, DummyLimiter())

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

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
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": "10.1000/xyz",
            }
            for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [{"scholar.PMID": pmid} for pmid in pmids]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", lambda *_, **__: {})

    gdd.fetch_pubmed_records(
        ["1"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert limiters["documents_openalex"].acquisitions == 1
    assert limiters["documents_crossref"].acquisitions == 1


def test_crossref_jobs_skip_missing_doi(monkeypatch: pytest.MonkeyPatch) -> None:
    """CrossRef jobs ignore records without a DOI when scheduling fetches."""

    class TrackingLimiter:
        def __init__(self) -> None:
            self.acquisitions = 0
            self._lock = threading.Lock()

        def acquire(self) -> None:
            with self._lock:
                self.acquisitions += 1

    limiters = {
        "documents_global": DummyLimiter(),
        "documents_crossref": TrackingLimiter(),
    }

    def fake_get_limiter(name: str, *_, **__) -> Any:
        return limiters.get(name, DummyLimiter())

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

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
        retry_cfg: Any | None = None,
        client: Any | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": "10.1000/xyz" if pmid != "1" else "",
            }
            for pmid in batch
        ]

    def fake_semantic_batch(
        session: Any,
        pmids: list[str],
        sleep: float,
        cfg: Any | None = None,
    ) -> list[dict[str, str]]:
        return [{"scholar.PMID": pmid} for pmid in pmids]

    seen_crossref: list[str] = []

    def fake_crossref(session: Any, doi: str, cfg: Any, limiter: Any) -> dict[str, str]:
        seen_crossref.append(doi)
        return {}

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_openalex", lambda *_, **__: {})
    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)

    gdd.fetch_pubmed_records(
        ["1", "2"],
        sleep=0.0,
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert limiters["documents_crossref"].acquisitions == 1
    assert seen_crossref == ["10.1000/xyz"]


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
        retry_cfg: Any | None = None,
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
        retry_cfg: Any | None = None,
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
        retry_cfg: Any | None = None,
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
        retry_cfg: Any | None = None,
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
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    document_export_postprocess_stub,
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
    monkeypatch.setattr(gdd, "build_quality_report", lambda *_, **__: {})
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(
        gdd, "analyze_table_quality", lambda df, table_name, **_: None
    )
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
    assert document_export_postprocess_stub


def test_finalise_export_accepts_series_chunk_size(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    document_export_postprocess_stub,
) -> None:
    """Chunk size provided as a scalar ``Series`` should be coerced to ``int``."""

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
        col_order: list[str] | None = None,
        chunksize: int,
        merge_chunksize: int,
        sort_chunksize: int,
        sep: str,
        encoding: str,
        **_: Any,
    ) -> Path:
        list(chunks)
        path.write_text("data")
        captured["chunksize"] = chunksize
        captured["merge_chunksize"] = merge_chunksize
        captured["sort_chunksize"] = sort_chunksize
        captured["key_cols"] = list(key_cols or [])
        captured["col_order"] = list(col_order or [])
        captured["sep"] = sep
        captured["encoding"] = encoding
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "build_quality_report", lambda *_, **__: {})
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda profiler, table_name: None)
    monkeypatch.setattr(gdd.DocumentsSchema, "validate", lambda frame, lazy=True: frame)

    exit_code = gdd._finalise_export(
        df,
        output,
        cfg,
        input_csv=tmp_path / "input.csv",
        key_columns=None,
        chunk_size=pd.Series([7]),
    )

    assert exit_code == 0
    assert captured["chunksize"] == 7
    assert captured["merge_chunksize"] == 7
    assert captured["sort_chunksize"] == 7
    assert captured["key_cols"] == ["ChEMBL.document_chembl_id"]
    assert document_export_postprocess_stub


def test_finalise_export_series_chunk_size_fallback(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    document_export_postprocess_stub,
) -> None:
    """Non-scalar ``Series`` chunk sizes should fall back to the default."""

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
        col_order: list[str] | None = None,
        chunksize: int,
        merge_chunksize: int,
        sort_chunksize: int,
        sep: str,
        encoding: str,
        **_: Any,
    ) -> Path:
        list(chunks)
        path.write_text("data")
        captured["chunksize"] = chunksize
        captured["merge_chunksize"] = merge_chunksize
        captured["sort_chunksize"] = sort_chunksize
        return path

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "build_quality_report", lambda *_, **__: {})
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda profiler, table_name: None)
    monkeypatch.setattr(gdd.DocumentsSchema, "validate", lambda frame, lazy=True: frame)

    exit_code = gdd._finalise_export(
        df,
        output,
        cfg,
        input_csv=tmp_path / "input.csv",
        key_columns=None,
        chunk_size=pd.Series([5, 6]),
    )

    assert exit_code == 0
    assert captured["chunksize"] == gdd._EXPORT_STREAM_CHUNK_SIZE
    assert captured["merge_chunksize"] == gdd._EXPORT_STREAM_CHUNK_SIZE
    assert captured["sort_chunksize"] == gdd._EXPORT_STREAM_CHUNK_SIZE
    assert document_export_postprocess_stub


def test_finalise_export_accepts_generator(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    document_export_postprocess_stub,
) -> None:
    """Generator inputs should be validated and written in order."""

    cfg = Config()
    frames = [
        pd.DataFrame(
            {
                "document_chembl_id": ["CHEMBL1"],
                "PubMed.PMID": ["101"],
            }
        ),
        pd.DataFrame(
            {
                "document_chembl_id": ["CHEMBL2"],
                "PubMed.PMID": ["102"],
            }
        ),
    ]
    output = tmp_path / "documents.csv"

    captured: dict[str, Any] = {"chunks": [], "validated": []}

    def fake_write_csv_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        *,
        cfg: Any,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **_: Any,
    ) -> Path:
        captured["chunks"] = [chunk.copy() for chunk in chunks]
        captured["key_cols"] = list(key_cols or [])
        return path

    def fake_validate(frame: pd.DataFrame, lazy: bool = True) -> pd.DataFrame:
        captured["validated"].append(frame.copy())
        return frame

    monkeypatch.setattr(gdd, "write_csv_chunks_deterministic", fake_write_csv_chunks)
    monkeypatch.setattr(gdd.DocumentsSchema, "validate", fake_validate)
    monkeypatch.setattr(gdd, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    def fake_build_report(data: Any) -> dict[str, Any]:
        rows = getattr(data, "rows_total", len(data)) if hasattr(data, "__len__") else getattr(data, "rows_total", 0)
        return {"rows": rows}

    monkeypatch.setattr(gdd, "build_quality_report", fake_build_report)
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)
    monkeypatch.setattr(
        gdd, "analyze_table_quality", lambda df, table_name, **_: None
    )

    exit_code = gdd._finalise_export(
        iter(frames),
        output,
        cfg,
        input_csv=tmp_path / "input.csv",
        key_columns=["document_chembl_id"],
    )

    assert exit_code == 0
    assert [chunk["PubMed.PMID"].tolist() for chunk in captured["chunks"]] == [
        ["101"],
        ["102"],
    ]
    assert [frame["PubMed.PMID"].tolist() for frame in captured["validated"]] == [
        ["101"],
        ["102"],
    ]
    assert document_export_postprocess_stub


def test_finalise_export_streaming_is_deterministic(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    document_export_postprocess_stub,
) -> None:
    """Streaming export should preserve deterministic ordering and payload."""

    cfg = Config()
    output = tmp_path / "documents.csv"

    def frame_generator() -> Iterator[pd.DataFrame]:
        yield pd.DataFrame(
            {
                "document_chembl_id": ["CHEMBL2"],
                "PubMed.PMID": ["102"],
            }
        )
        yield pd.DataFrame(
            {
                "document_chembl_id": ["CHEMBL1"],
                "PubMed.PMID": ["101"],
            }
        )

    captured: dict[str, Any] = {}

    monkeypatch.setattr(gdd.DocumentsSchema, "validate", lambda frame, lazy=True: frame)
    monkeypatch.setattr(gdd, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)

    def fake_build_quality_report(data: Any) -> dict[str, Any]:
        captured["quality_input"] = data
        if isinstance(data, pd.DataFrame):
            captured["quality_df"] = data.copy()
            return {"rows": len(data)}
        rows = getattr(data, "rows_total", 0)
        return {"rows": rows}

    monkeypatch.setattr(gdd, "build_quality_report", fake_build_quality_report)
    monkeypatch.setattr(gdd, "save_quality_report", lambda report, path: path)

    def fake_analyze_table_quality(data: Any, table_name: str) -> None:
        captured["quality_table_name"] = table_name
        captured["quality_analyzer_input"] = data
        return None

    monkeypatch.setattr(gdd, "analyze_table_quality", fake_analyze_table_quality)

    exit_code = gdd._finalise_export(
        frame_generator(),
        output,
        cfg,
        input_csv=tmp_path / "input.csv",
        key_columns=["PubMed.PMID"],
        chunk_size=1,
    )

    exported = pd.read_csv(output, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding)

    assert exit_code == 0
    assert exported["PubMed.PMID"].astype(str).tolist() == ["101", "102"]
    quality_input = captured["quality_input"]
    assert getattr(quality_input, "rows_total", 0) == 2
    analyzer_input = captured["quality_analyzer_input"]
    quality_report, _ = analyzer_input.build(captured["quality_table_name"])
    assert quality_report["PubMed.PMID"].astype(str).tolist() == ["101", "102"]
    assert document_export_postprocess_stub


@pytest.mark.parametrize("context_position", ["suffix", "prefix"])
def test_fetch_pubmed_records_accepts_executor_context(
    monkeypatch: pytest.MonkeyPatch,
    context_position: str,
) -> None:
    """Executor passing an internal context argument should be ignored."""

    pmids = ["1", "2"]

    class DummyLimiter:
        def __init__(self) -> None:
            self.burst = 1

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


def test_fetch_pubmed_records_parallel_enrichment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OpenAlex and CrossRef lookups should run in parallel and preserve order."""

    pmids = ["100", "200", "300"]

    monkeypatch.setattr(rl, "sleep", lambda _delay: None)

    def fake_pubmed_batch(
        _session: Any,
        batch: Sequence[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"10.1234/{pmid}",
                "PubMed.ArticleTitle": f"Article {pmid}",
                "PubMed.PublicationType": "",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)

    def fake_semantic_batch(
        _session: Any,
        batch: Sequence[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1234/{pmid}",
                "scholar.PublicationTypes": "",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_args, **_kwargs: {
            "scholar.PMID": "",
            "scholar.DOI": "",
            "scholar.PublicationTypes": "",
            "scholar.Error": "",
        },
    )

    openalex_threads: list[str] = []
    openalex_completion: list[str] = []

    def fake_openalex(
        _session: Any,
        pmid: str,
        cfg: OpenAlexCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        openalex_threads.append(threading.current_thread().name)
        if pmid == "100":
            time.sleep(0.01)
        openalex_completion.append(pmid)
        return {
            "OpenAlex.PublicationTypes": "",
            "OpenAlex.TypeCrossref": "",
            "OpenAlex.Genre": "",
            "OpenAlex.Id": pmid,
            "OpenAlex.Venue": f"Venue {pmid}",
            "OpenAlex.MeshDescriptors": "",
            "OpenAlex.MeshQualifiers": "",
            "OpenAlex.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)

    crossref_threads: list[str] = []
    crossref_completion: list[str] = []

    def fake_crossref(
        _session: Any,
        doi: str,
        cfg: CrossRefCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        crossref_threads.append(threading.current_thread().name)
        pmid = doi.rsplit("/", 1)[-1] if doi else ""
        if pmid == "200":
            time.sleep(0.01)
        crossref_completion.append(pmid)
        return {
            "crossref.Type": "",
            "crossref.Subtype": "",
            "crossref.Title": f"Title {pmid}",
            "crossref.Subtitle": "",
            "crossref.Subject": "",
            "crossref.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=3,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert df["OpenAlex.Venue"].tolist() == [f"Venue {pmid}" for pmid in pmids]
    assert df["crossref.Title"].tolist() == [f"Title {pmid}" for pmid in pmids]

    assert openalex_completion != pmids
    assert crossref_completion != pmids

    assert openalex_threads
    assert crossref_threads
    assert all(name.startswith("ThreadPoolExecutor") for name in openalex_threads)
    assert all(name.startswith("ThreadPoolExecutor") for name in crossref_threads)


def test_fetch_pubmed_records_reuses_service_sessions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Service lookups should reuse sessions within a batch."""

    pmids = ["100", "200", "300"]
    session_labels: list[str] = []

    class DummySession:
        def __init__(self, label: str) -> None:
            self.label = label

        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

    def fake_session_with_retry(*_args: object, **_kwargs: object) -> DummySession:
        label = f"session-{len(session_labels)}"
        session_labels.append(label)
        return DummySession(label)

    monkeypatch.setattr(gdd, "session_with_retry", fake_session_with_retry)

    openalex_session_labels: list[str] = []
    crossref_session_labels: list[str] = []

    def fake_openalex_session(*_args: object, **_kwargs: object) -> DummySession:
        label = f"openalex-{len(openalex_session_labels)}"
        openalex_session_labels.append(label)
        return DummySession(label)

    def fake_crossref_session(*_args: object, **_kwargs: object) -> DummySession:
        label = f"crossref-{len(crossref_session_labels)}"
        crossref_session_labels.append(label)
        return DummySession(label)

    monkeypatch.setattr(gdd, "openalex_session", fake_openalex_session)
    monkeypatch.setattr(gdd, "crossref_session", fake_crossref_session)

    def fake_pubmed_batch(
        session: DummySession,
        batch: Sequence[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        assert session.label == "session-0"
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"10.1234/{pmid}",
                "PubMed.ArticleTitle": f"Article {pmid}",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)

    def fake_semantic_batch(
        session: DummySession,
        batch: Sequence[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        assert session.label == "session-0"
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1234/{pmid}",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)

    openalex_sessions: list[str] = []

    def fake_openalex(
        session: DummySession,
        pmid: str,
        cfg: OpenAlexCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        openalex_sessions.append(session.label)
        return {
            "OpenAlex.Id": pmid,
            "OpenAlex.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_openalex)

    crossref_sessions: list[str] = []

    def fake_crossref(
        session: DummySession,
        doi: str,
        cfg: CrossRefCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        crossref_sessions.append(session.label)
        return {
            "crossref.Title": doi,
            "crossref.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_crossref)
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    class ImmediateFuture:
        def __init__(self, value: list[dict[str, str]]) -> None:
            self._value = value

        def result(self) -> list[dict[str, str]]:  # pragma: no cover - trivial
            return self._value

    class ImmediateExecutor:
        def __init__(self, *args: object, **kwargs: object) -> None:
            return None

        def __enter__(self) -> "ImmediateExecutor":  # pragma: no cover - trivial
            return self

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            return None

        def submit(
            self,
            fn: Callable[[Sequence[str]], list[dict[str, str]]],
            batch: Sequence[str],
        ) -> ImmediateFuture:
            return ImmediateFuture(fn(batch))

    monkeypatch.setattr(gdd, "ThreadPoolExecutor", ImmediateExecutor)
    monkeypatch.setattr(gdd, "as_completed", lambda futures: list(futures))

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(burst=1, rps=1),
        crossref_cfg=CrossRefCfg(burst=1, rps=1),
        max_workers=1,
        batch_size=3,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert len(session_labels) == 1
    assert len(openalex_session_labels) == 1
    assert len(crossref_session_labels) == 1
    assert openalex_sessions and len(openalex_sessions) == len(pmids)
    assert openalex_sessions == [openalex_sessions[0]] * len(openalex_sessions)
    assert crossref_sessions and len(crossref_sessions) == len(pmids)
    assert crossref_sessions == [crossref_sessions[0]] * len(crossref_sessions)


@responses.activate
def test_fetch_pubmed_records_reuses_sessions_across_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Worker threads should reuse HTTP sessions across multiple batches."""

    pmids = [str(100 + i) for i in range(6)]

    for pmid in pmids:
        responses.add(
            responses.GET,
            f"https://pubmed.test/{pmid}",
            json={"pmid": pmid},
            status=200,
        )
        responses.add(
            responses.GET,
            f"https://openalex.test/{pmid}",
            json={"id": pmid},
            status=200,
        )
        responses.add(
            responses.GET,
            f"https://crossref.test/10.1234/{pmid}",
            json={"doi": pmid},
            status=200,
        )

    class RecordingContext:
        def __init__(self, label: str) -> None:
            self.label = label
            self.session = requests.Session()
            self.closed = False
            self.session_id = id(self.session)

        def __enter__(self) -> requests.Session:  # pragma: no cover - trivial
            return self.session

        def __exit__(self, *exc: object) -> None:  # pragma: no cover - trivial
            self.closed = True
            self.session.close()

    contexts: list[RecordingContext] = []

    def _make_context(label: str) -> RecordingContext:
        ctx = RecordingContext(label)
        contexts.append(ctx)
        return ctx

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: _make_context("pubmed"))
    monkeypatch.setattr(gdd, "openalex_session", lambda *_, **__: _make_context("openalex"))
    monkeypatch.setattr(gdd, "crossref_session", lambda *_, **__: _make_context("crossref"))
    monkeypatch.setattr(gdd, "get_limiter", lambda *_, **__: DummyLimiter())

    pubmed_session_ids: set[int] = set()
    openalex_session_ids: set[int] = set()
    crossref_session_ids: set[int] = set()

    def fake_pubmed_batch(
        session: requests.Session,
        batch: Sequence[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        records: list[dict[str, str]] = []
        for pmid in batch:
            pubmed_session_ids.add(id(session))
            session.get(
                f"https://pubmed.test/{pmid}",
                headers={"X-Session": str(id(session))},
                timeout=0.1,
            )
            records.append(
                {
                    "PubMed.PMID": pmid,
                    "PubMed.DOI": f"10.1234/{pmid}",
                    "PubMed.ArticleTitle": f"Article {pmid}",
                }
            )
        return records

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)

    def fake_semantic_batch(
        _session: requests.Session,
        batch: Sequence[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1234/{pmid}",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_args, **_kwargs: {
            "scholar.PMID": "",
            "scholar.DOI": "",
            "scholar.Error": "unreachable",
        },
    )

    def fake_fetch_openalex(
        session: requests.Session,
        pmid: str,
        cfg: OpenAlexCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        openalex_session_ids.add(id(session))
        session.get(
            f"https://openalex.test/{pmid}",
            headers={"X-Session": str(id(session))},
            timeout=0.1,
        )
        return {
            "OpenAlex.PMID": pmid,
            "OpenAlex.DOI": f"10.1234/{pmid}",
            "OpenAlex.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_openalex", fake_fetch_openalex)

    def fake_fetch_crossref(
        session: requests.Session,
        doi: str,
        cfg: CrossRefCfg,
        limiter: rl.RateLimiter | None = None,
    ) -> dict[str, str]:
        crossref_session_ids.add(id(session))
        session.get(
            f"https://crossref.test/{doi}",
            headers={"X-Session": str(id(session))},
            timeout=0.1,
        )
        return {
            "crossref.DOI": doi,
            "crossref.Title": f"Title {doi}",
            "crossref.Error": "",
        }

    monkeypatch.setattr(gdd.ocl, "fetch_crossref", fake_fetch_crossref)

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=2,
        batch_size=1,
    )

    assert df["PubMed.PMID"].tolist() == pmids

    assert len(pubmed_session_ids) <= 2
    assert len(openalex_session_ids) <= 2
    assert len(crossref_session_ids) <= 2

    pubmed_headers = {
        call.request.headers.get("X-Session")
        for call in responses.calls
        if call.request.url.startswith("https://pubmed.test/")
    }
    openalex_headers = {
        call.request.headers.get("X-Session")
        for call in responses.calls
        if call.request.url.startswith("https://openalex.test/")
    }
    crossref_headers = {
        call.request.headers.get("X-Session")
        for call in responses.calls
        if call.request.url.startswith("https://crossref.test/")
    }

    assert len(pubmed_headers) <= 2
    assert len(openalex_headers) <= 2
    assert len(crossref_headers) <= 2

    assert all(ctx.closed for ctx in contexts)


def test_fetch_pubmed_records_reuses_service_executors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OpenAlex and CrossRef executors should be shared across batches."""

    pmids = ["100", "200", "300", "400"]
    creations: list[int | None] = []

    orig_executor = gdd.ThreadPoolExecutor

    class TrackingExecutor(orig_executor):
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            max_workers = kwargs.get("max_workers")
            if max_workers is None and args:
                max_workers = args[0]
            creations.append(max_workers)
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(gdd, "ThreadPoolExecutor", TrackingExecutor)
    monkeypatch.setattr(
        gdd, "get_limiter", lambda _name, _rps, burst=None: DummyLimiter(burst or 1)
    )

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *_exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_, **__: DummySession())

    def fake_pubmed_batch(
        _session: Any,
        batch: Sequence[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": f"10.1234/{pmid}",
                "PubMed.ArticleTitle": f"Article {pmid}",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)

    def fake_semantic_batch(
        _session: Any,
        batch: Sequence[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": f"10.1234/{pmid}",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(gdd.ssl, "fetch_semantic_scholar_batch", fake_semantic_batch)

    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_args, **_kwargs: {"scholar.PMID": "", "scholar.DOI": ""},
    )

    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda *_args, **_kwargs: {
            "OpenAlex.Venue": "Venue",
            "OpenAlex.Error": "",
        },
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda *_args, **_kwargs: {
            "crossref.Title": "Title",
            "crossref.Error": "",
        },
    )

    openalex_cfg = OpenAlexCfg()
    crossref_cfg = CrossRefCfg()
    max_workers = 1

    df = gdd.fetch_pubmed_records(
        pmids,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=openalex_cfg,
        crossref_cfg=crossref_cfg,
        max_workers=max_workers,
        batch_size=2,
    )

    assert df["PubMed.PMID"].tolist() == pmids
    assert len(creations) == 3

    expected_downstream_capacity = max(
        1,
        min(max_workers, openalex_cfg.burst, crossref_cfg.burst),
    )
    assert creations[0] == expected_downstream_capacity


def test_fetch_pubmed_records_generator_batches(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Generator mode should yield DataFrame batches in submission order."""

    pmids = ["100", "200", "300", "400"]
    cfg = Config()

    def fake_get_limiter(*_: object, **__: object) -> DummyLimiter:
        return DummyLimiter(burst=10)

    monkeypatch.setattr(gdd, "get_limiter", fake_get_limiter)

    def fake_as_completed(futures: Iterable[Future]) -> Iterator[Future]:
        snapshot = list(futures)
        for future in reversed(snapshot):
            yield future

    monkeypatch.setattr(gdd, "as_completed", fake_as_completed)

    class ImmediateExecutor:
        def __init__(self, *_: object, **__: object) -> None:
            return None

        def __enter__(self) -> ImmediateExecutor:
            return self

        def __exit__(self, *exc: object) -> None:
            return None

        def submit(self, func: Any, *args: object, **__: object) -> Future:
            _ = func
            batch = list(args[0]) if args else []
            future: Future = Future()
            records = [{"PubMed.PMID": pmid} for pmid in batch]
            future.set_result(records)
            return future

    monkeypatch.setattr(gdd, "ThreadPoolExecutor", ImmediateExecutor)

    generator = gdd.fetch_pubmed_records(
        pmids,
        cfg,
        return_generator=True,
        batch_size=1,
        max_workers=2,
    )

    assert isinstance(generator, Iterator)
    frames = list(generator)
    assert [frame["PubMed.PMID"].tolist() for frame in frames] == [
        ["100"],
        ["200"],
        ["300"],
        ["400"],
    ]


def test_fetch_pubmed_records_uses_pending_helper_minimally(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure ``_iter_pending`` only triggers ``as_completed`` when required."""

    pmids = ["100", "200"]
    cfg = Config()

    monkeypatch.setattr(
        gdd, "get_limiter", lambda *_args, **_kwargs: DummyLimiter(burst=1)
    )

    class ImmediateExecutor:
        def __init__(self, *_: object, **__: object) -> None:
            return None

        def __enter__(self) -> ImmediateExecutor:
            return self

        def __exit__(self, *_exc: object) -> None:
            return None

        def submit(
            self, func: Callable[..., list[dict[str, str]]], *args: object, **kwargs: object
        ) -> Future:
            future: Future = Future()
            try:
                result = func(*args, **kwargs)
            except Exception as exc:  # pragma: no cover - defensive
                future.set_exception(exc)
            else:
                future.set_result(result)
            return future

    monkeypatch.setattr(gdd, "ThreadPoolExecutor", ImmediateExecutor)

    class DummySession:
        def __enter__(self) -> DummySession:  # pragma: no cover - trivial
            return self

        def __exit__(self, *_exc: object) -> None:  # pragma: no cover - trivial
            return None

    monkeypatch.setattr(gdd, "session_with_retry", lambda *_args, **_kwargs: DummySession())

    def fake_pubmed_batch(
        _session: object,
        batch: Sequence[str],
        _sleep: float,
        cfg: PubMedCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "PubMed.PMID": pmid,
                "PubMed.DOI": "",
                "PubMed.ArticleTitle": f"Article {pmid}",
            }

            for pmid in batch
        ]

    monkeypatch.setattr(gdd.pl, "fetch_pubmed_batch", fake_pubmed_batch)

    def fake_semantic_batch(

        _session: object,
        batch: Sequence[str],
        _sleep: float,
        cfg: SemanticScholarCfg | None = None,
    ) -> list[dict[str, str]]:
        return [
            {
                "scholar.PMID": pmid,
                "scholar.DOI": "",
                "scholar.Error": "",
            }
            for pmid in batch
        ]

    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar_batch",
        fake_semantic_batch,
    )
    monkeypatch.setattr(
        gdd.ssl,
        "fetch_semantic_scholar",
        lambda *_args, **_kwargs: {
            "scholar.PMID": "",
            "scholar.DOI": "",
            "scholar.Error": "",
        },
    )

    monkeypatch.setattr(
        gdd.ocl,
        "fetch_openalex",
        lambda *_args, **_kwargs: {"OpenAlex.Error": ""},
    )
    monkeypatch.setattr(
        gdd.ocl,
        "fetch_crossref",
        lambda *_args, **_kwargs: {"crossref.Error": ""},
    )

    call_count = 0
    orig_as_completed = gdd.as_completed

    def spy_as_completed(futures: Iterable[Future]) -> Iterator[Future]:
        nonlocal call_count
        call_count += 1
        return orig_as_completed(futures)

    monkeypatch.setattr(gdd, "as_completed", spy_as_completed)

    df = gdd.fetch_pubmed_records(
        pmids,
        cfg,
        sleep=0.0,
        pubmed_cfg=PubMedCfg(),
        semantic_scholar_cfg=SemanticScholarCfg(),
        openalex_cfg=OpenAlexCfg(),
        crossref_cfg=CrossRefCfg(),
        max_workers=1,
        batch_size=1,
    )

    assert call_count == 2
    assert df["PubMed.PMID"].tolist() == pmids


def test_prepare_export_frame_merges_prefixed_columns() -> None:
    """Existing prefixed columns should absorb their unprefixed counterparts."""

    frame = pd.DataFrame(
        {
            "ChEMBL.title": ["kept", ""],
            "title": ["ignored", "fallback"],
            "ChEMBL.abstract": ["", "pref"],
            "abstract": ["primary", ""],
        }
    )

    result = gdd._prepare_export_frame(frame)

    assert not result.columns.duplicated().any()
    assert result.loc[0, "ChEMBL.title"] == "kept"
    assert result.loc[1, "ChEMBL.title"] == "fallback"
    assert result.loc[0, "ChEMBL.abstract"] == "primary"
    assert result.loc[1, "ChEMBL.abstract"] == "pref"

