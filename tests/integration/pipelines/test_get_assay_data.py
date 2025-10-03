from __future__ import annotations

import argparse
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd
import pytest

import library.cli_utils as cli_utils
from library.integration import chembl_library as cl
from library import io
from library.config import Config
from scripts import get_assay_data as gas


def test_run_chembl_orders_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Output columns should follow `AssaysSchema` order with extras alphabetic."""
    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nA1\n")
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["A1"]))

    df = pd.DataFrame(
        {
            "zzz": ["z"],
            "document_chembl_id": ["D1"],
            "aaa": ["a"],
            "assay_chembl_id": ["A1"],
            "target_chembl_id": ["T1"],
        }
    )
    monkeypatch.setattr(cl, "get_assays", lambda *_, **__: df)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda df: df)
    monkeypatch.setattr(
        gas, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    rc = gas.run_chembl(cfg, args)
    assert rc == 0
    out_df = pd.read_csv(args.output_csv)
    assert list(out_df.columns) == [
        "assay_chembl_id",
        "document_chembl_id",
        "target_chembl_id",
        "pipeline_version",
        "timestamp_utc",
        "aaa",
        "zzz",
    ]


def test_run_chembl_streams_chunks_and_applies_hooks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Chunk iterator remains lazy and includes metadata/validator output."""

    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nA1\nA2\n")
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["A1", "A2"]))

    df = pd.DataFrame(
        {
            "assay_chembl_id": ["A1", "A2"],
            "document_chembl_id": ["D1", "D2"],
        }
    )

    monkeypatch.setattr(cl, "get_assays", lambda *_, **__: df)

    hook_calls: dict[str, int] = {"postprocess": 0, "normalize": 0, "metadata": 0}

    def fake_postprocess(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["postprocess"] += 1
        return frame.assign(postprocessed=True)

    def fake_normalize(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["normalize"] += 1
        assert "postprocessed" in frame.columns
        return frame.assign(normalized=True)

    def fake_metadata(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["metadata"] += 1
        assert "normalized" in frame.columns
        return frame.assign(pipeline_version="1.0", timestamp_utc="now")

    monkeypatch.setattr(gas.ap, "postprocess_assays", fake_postprocess)
    monkeypatch.setattr(gas, "normalize_assays", fake_normalize)
    monkeypatch.setattr(gas, "add_pipeline_metadata", fake_metadata)

    class _Result:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    validator_calls = {"count": 0}

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> _Result:
        validator_calls["count"] += 1
        assert "timestamp_utc" in frame.columns
        return _Result(frame)

    monkeypatch.setattr(gas, "validate_assays", fake_validate)
    monkeypatch.setattr(
        gas, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, object] = {}

    def fake_write(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        *,
        key_cols: Sequence[str],
        col_order: Sequence[str],
        chunksize: int,
        sort_chunksize: int,
        sep: str,
        encoding: str,
        cfg: Config,
        **kwargs: object,
    ) -> Path:
        assert not isinstance(chunks, list)
        frames = list(chunks)
        captured["chunk_count"] = len(frames)
        captured["columns"] = [set(frame.columns) for frame in frames]
        captured["key_cols"] = list(key_cols)
        return destination

    monkeypatch.setattr(gas, "write_csv_chunks_deterministic", fake_write)

    rc = gas.run_chembl(cfg, args)
    assert rc == 0
    assert hook_calls == {"postprocess": 1, "normalize": 1, "metadata": 1}
    assert validator_calls["count"] == 1
    assert captured["chunk_count"] == 1
    assert all(
        {"postprocessed", "normalized", "pipeline_version", "timestamp_utc"}.issubset(
            cols
        )
        for cols in captured["columns"]
    )
    assert captured["key_cols"] == ["assay_chembl_id"]


def test_run_chembl_processes_multiple_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Assay chunks are fetched lazily and processed individually."""

    cfg.sources.chembl.pipelines.assay.batch_size = 1
    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nA1\nA2\nA3\n")
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    ids = ["A1", "A2", "A3"]
    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(ids))

    captured_calls: list[list[str]] = []

    def fake_get_assays(chunk_ids: Iterable[str], **_: object) -> pd.DataFrame:
        ids_list = list(chunk_ids)
        assert len(ids_list) == 1
        captured_calls.append(ids_list)
        assay_id = ids_list[0]
        return pd.DataFrame(
            {
                "assay_chembl_id": [assay_id],
                "document_chembl_id": [f"D{assay_id}"],
            }
        )

    monkeypatch.setattr(cl, "get_assays", fake_get_assays)

    hook_calls: dict[str, int] = {"postprocess": 0, "normalize": 0, "metadata": 0}

    def fake_postprocess(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["postprocess"] += 1
        return frame.assign(postprocessed=True)

    def fake_normalize(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["normalize"] += 1
        assert "postprocessed" in frame.columns
        return frame.assign(normalized=True)

    def fake_metadata(frame: pd.DataFrame) -> pd.DataFrame:
        hook_calls["metadata"] += 1
        assert "normalized" in frame.columns
        return frame.assign(pipeline_version="1.0", timestamp_utc="now")

    monkeypatch.setattr(gas.ap, "postprocess_assays", fake_postprocess)
    monkeypatch.setattr(gas, "normalize_assays", fake_normalize)
    monkeypatch.setattr(gas, "add_pipeline_metadata", fake_metadata)

    class _Result:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    validator_calls = 0

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> _Result:
        nonlocal validator_calls
        validator_calls += 1
        assert "timestamp_utc" in frame.columns
        assert return_result is True
        return _Result(frame)

    monkeypatch.setattr(gas, "validate_assays", fake_validate)
    monkeypatch.setattr(
        gas, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, object] = {}

    def fake_write(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        *,
        key_cols: Sequence[str],
        col_order: Sequence[str],
        chunksize: int,
        sort_chunksize: int,
        sep: str,
        encoding: str,
        cfg: Config,
        **kwargs: object,
    ) -> Path:
        assert not isinstance(chunks, list)
        frames = list(chunks)
        captured["chunk_count"] = len(frames)
        captured["rows"] = [frame["assay_chembl_id"].tolist() for frame in frames]
        captured["key_cols"] = list(key_cols)
        return destination

    monkeypatch.setattr(gas, "write_csv_chunks_deterministic", fake_write)

    rc = gas.run_chembl(cfg, args)

    assert rc == 0
    assert captured_calls == [["A1"], ["A2"], ["A3"]]
    assert hook_calls == {"postprocess": 3, "normalize": 3, "metadata": 3}
    assert validator_calls == 3
    assert captured["chunk_count"] == 3
    assert captured["rows"] == [["A1"], ["A2"], ["A3"]]
    assert captured["key_cols"] == ["assay_chembl_id"]


def test_run_chembl_sorts_output_deterministically(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Rows are sorted by key columns when writing output."""

    input_csv = tmp_path / "assays.csv"
    input_csv.write_text("assay_chembl_id\nA2\nA1\n")
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["A2", "A1"]))

    df = pd.DataFrame(
        {
            "assay_chembl_id": ["A2", "A1"],
            "document_chembl_id": ["D2", "D1"],
        }
    )

    monkeypatch.setattr(cl, "get_assays", lambda *_, **__: df)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda frame: frame)
    monkeypatch.setattr(
        gas, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    rc = gas.run_chembl(cfg, args)
    assert rc == 0

    result = pd.read_csv(args.output_csv)
    assert list(result["assay_chembl_id"]) == ["A1", "A2"]
