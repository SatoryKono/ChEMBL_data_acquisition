import argparse
from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from library.pipelines.document import runner


class _StubPipeline:
    def __init__(self) -> None:
        self.limit_calls: list[tuple[list[str], int]] = []
        self.fetch_calls: list[tuple[list[str], dict[str, object]]] = []
        self.read_calls: list[tuple[Path, int]] = []

    def build_missing_input_context(self, path: Path) -> dict[str, str]:
        return {"path": str(path)}

    def limit_iterable(self, iterable, limit: int):
        values = list(iterable)
        limited = values[:limit]
        self.limit_calls.append((values, limit))

        def _counter() -> int:
            return len(limited)

        return iter(limited), _counter

    def fetch_pubmed_records(self, pmids, **kwargs):
        pmid_list = list(pmids)
        self.fetch_calls.append((pmid_list, dict(kwargs)))
        rows = [
            {
                "document_chembl_id": f"CHEMBL{pmid}",
                "PubMed.PMID": pmid,
                "value": f"meta-{pmid}",
            }
            for pmid in pmid_list
        ]
        yield pd.DataFrame(rows)

    def read_csv_chunks(self, path: Path, *, cfg: Config, chunk_size: int):
        self.read_calls.append((Path(path), chunk_size))
        frame = pd.DataFrame(
            {
                "document_chembl_id": ["CHEMBL101"],
                "PubMed.PMID": ["101"],
                "merged": ["yes"],
            }
        )
        return iter([frame])

    def apply_fallback_doi(self, frame: pd.DataFrame, **kwargs) -> pd.DataFrame:
        return frame

    def build_fallback_doi_map(self, frame: pd.DataFrame, **kwargs) -> dict[str, str]:
        return {}


@pytest.fixture()
def pipeline_stub() -> _StubPipeline:
    return _StubPipeline()


def test_run_pubmed__streams_frames(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pipeline_stub: _StubPipeline
) -> None:
    pmids = ["100", "101", "102"]
    monkeypatch.setattr(
        runner.io,
        "read_ids",
        lambda path, column, cfg: iter(pmids),
    )

    captured: dict[str, object] = {}

    def _fake_finalise(df, output, cfg, *, input_csv, key_columns, chunk_size):
        frames = list(df)
        captured["frames"] = frames
        captured["output"] = output
        captured["key_columns"] = key_columns
        captured["chunk_size"] = chunk_size
        return 0

    monkeypatch.setattr(runner, "_finalise_export", _fake_finalise)
    monkeypatch.setattr(runner, "normalize_documents", lambda frame: frame)

    output_csv = tmp_path / "pubmed.csv"
    args = argparse.Namespace(
        input_csv=tmp_path / "pmids.csv",
        final_out=output_csv,
        output_csv=output_csv,
        limit=2,
        offset=1,
        batch_size=1,
        sleep=0.0,
        workers=1,
        fallback_doi_enabled=False,
        fallback_doi_path=None,
        fallback_doi_overwrite=False,
    )

    exit_code = runner.run_pubmed(cfg, args, pipeline=pipeline_stub)

    assert exit_code == 0
    assert captured["output"] == output_csv
    frames = captured["frames"]
    assert len(frames) == 1
    assert list(frames[0]["PubMed.PMID"]) == ["101", "102"]
    assert pipeline_stub.limit_calls == [(["101", "102"], 2)]
    assert pipeline_stub.fetch_calls[0][0] == ["101", "102"]


def test_run_chembl__writes_finalised_data(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pipeline_stub: _StubPipeline
) -> None:
    chembl_ids = ["CHEMBL1", "CHEMBL2"]
    monkeypatch.setattr(
        runner.io,
        "read_ids",
        lambda path, column, cfg: iter(chembl_ids),
    )

    captured: dict[str, object] = {}

    def _fake_finalise(df, output, cfg, *, input_csv, key_columns, chunk_size):
        captured["df"] = df
        captured["output"] = output
        captured["chunk_size"] = chunk_size
        return 0

    monkeypatch.setattr(runner, "_finalise_export", _fake_finalise)
    monkeypatch.setattr(runner, "normalize_documents", lambda frame: frame)

    recorded_ids: list[list[str]] = []

    def _fake_get_documents(ids, **kwargs):
        collected = list(ids)
        recorded_ids.append(collected)
        return pd.DataFrame(
            {
                "document_chembl_id": collected,
                "doi": ["10.1000/example" for _ in collected],
            }
        )

    monkeypatch.setattr(runner.cl, "get_documents", _fake_get_documents)

    class _Client:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(runner, "ChemblClient", lambda *args, **kwargs: _Client())

    output_csv = tmp_path / "chembl.csv"
    args = argparse.Namespace(
        input_csv=tmp_path / "ids.csv",
        final_out=output_csv,
        output_csv=output_csv,
        limit=1,
        offset=0,
        chunk_size=1,
        timeout=1.0,
    )

    exit_code = runner.run_chembl(cfg, args, pipeline=pipeline_stub)

    assert exit_code == 0
    assert captured["output"] == output_csv
    assert recorded_ids == [["CHEMBL1"]]
    assert captured["chunk_size"] == 1


def test_run_all__merges_sources(
    cfg: Config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, pipeline_stub: _StubPipeline
) -> None:
    chembl_ids = ["CHEMBL1", "CHEMBL2"]
    monkeypatch.setattr(
        runner.io,
        "read_ids",
        lambda path, column, cfg: iter(chembl_ids),
    )

    def _fake_get_documents(ids, **kwargs):
        collected = list(ids)
        return pd.DataFrame(
            {
                "document_chembl_id": collected,
                "pubmed_id": [101, 0],
                "doi": ["10.1000/example", ""],
            }
        )

    monkeypatch.setattr(runner.cl, "get_documents", _fake_get_documents)

    def _fake_write(chunks, path, **kwargs):
        data_frames = list(chunks)
        if data_frames:
            combined = pd.concat(data_frames, ignore_index=True)
        else:
            combined = pd.DataFrame()
        combined.to_csv(path, index=False)
        return path

    monkeypatch.setattr(runner, "write_csv_chunks_deterministic", _fake_write)
    monkeypatch.setattr(runner, "merge_with_chembl", lambda docs, meta: docs.merge(next(meta, pd.DataFrame()), on="document_chembl_id", how="left"))
    monkeypatch.setattr(runner.dp, "postprocess_documents", lambda frame: frame)
    monkeypatch.setattr(runner, "normalize_documents", lambda frame: frame)

    captured: dict[str, object] = {}

    def _fake_finalise(df, output, cfg, *, input_csv, key_columns, chunk_size):
        captured["df"] = df
        captured["output"] = output
        captured["chunk_size"] = chunk_size
        return 0

    monkeypatch.setattr(runner, "_finalise_export", _fake_finalise)

    class _Client:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(runner, "ChemblClient", lambda *args, **kwargs: _Client())

    output_csv = tmp_path / "all.csv"
    args = argparse.Namespace(
        input_csv=tmp_path / "ids.csv",
        final_out=output_csv,
        output_csv=output_csv,
        limit=None,
        offset=0,
        chembl_chunk_size=2,
        chembl_timeout=1.0,
        pubmed_workers=1,
        pubmed_batch_size=2,
        pubmed_sleep=0.0,
        fallback_doi_enabled=False,
        fallback_doi_path=None,
        fallback_doi_overwrite=False,
    )

    exit_code = runner.run_all(cfg, args, pipeline=pipeline_stub)

    assert exit_code == 0
    assert captured["output"] == output_csv
    assert isinstance(captured["df"], pd.DataFrame)
    assert "merged" in captured["df"].columns
    assert 101 in captured["df"]["pubmed_id"].tolist()
    assert pipeline_stub.fetch_calls[0][0] == ["101", "0"]
    assert pipeline_stub.read_calls
