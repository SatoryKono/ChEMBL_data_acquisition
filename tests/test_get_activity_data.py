from __future__ import annotations

import argparse
import time
from collections.abc import Iterable
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import parse_qs, urlparse

import pandas as pd
import pytest

import library.cli_utils as cli_utils
from library import chembl_library as cl
from library import io
from library import rate_limiter as rl
from library.config import Config
from schemas import ActivitiesSchema
from scripts import get_activity_data as gad


def test_run_chembl_respects_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg.activity.limit = 2
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    read_ids: list[str] = []

    def fake_read_ids(*_, **__):
        for item in ["1", "2", "3"]:
            read_ids.append(item)
            yield item

    monkeypatch.setattr(io, "read_ids", fake_read_ids)

    captured: dict[str, object] = {}

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        data = list(ids)
        captured["ids"] = data
        captured["extra_columns"] = kwargs.get("extra_columns")
        return pd.DataFrame(
            {
                "activity_id": data,
                "assay_chembl_id": [f"A{i}" for i in data],
                "molecule_chembl_id": [f"M{i}" for i in data],
                "standard_value": [1.0] * len(data),
                "standard_type": ["IC50"] * len(data),
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)

    monkeypatch.setattr(
        gad,
        "write_csv_chunks_deterministic",
        lambda df, output, *, key_cols, col_order, chunksize, sort_chunksize, sep, encoding, cfg, **__: output,
    )
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    assert captured["ids"] == ["1", "2"]
    assert read_ids == ["1", "2"]
    assert captured["extra_columns"] == ["action_type"]


def test_run_chembl_limit_dry_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n2\n3\n")

    cfg.activity.limit = 2
    cfg.activity.dry_run = True
    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    called = {"read": False}

    def fake_read_ids(*_, **__):
        called["read"] = True
        return iter([])

    monkeypatch.setattr(io, "read_ids", fake_read_ids)

    rc = gad.run_chembl(cfg, args)
    assert rc == 0
    assert called["read"] is False


def test_run_chembl_column_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure schema columns precede alphabetically sorted extras."""
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["1"]))

    df = pd.DataFrame(
        [
            {
                "standard_value": 1.0,
                "assay_description": "desc",
                "activity_id": "1",
                "molecule_chembl_id": "CHEMBL1",
                "bao_format": "A",
                "standard_type": "IC50",
                "target_id": "T1",
                "assay_chembl_id": "A1",
            }
        ]
    )

    monkeypatch.setattr(cl, "get_activities", lambda *_, **__: df)
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, list[str]] = {}

    def fake_write_csv(
        df,
        output,
        *,
        key_cols,
        col_order,
        chunksize,
        sort_chunksize,
        sep,
        encoding,
        cfg,
        **__,
    ) -> Path:
        captured["col_order"] = list(col_order or [])
        return output

    monkeypatch.setattr(gad, "write_csv_chunks_deterministic", fake_write_csv)

    rc = gad.run_chembl(cfg, args)
    assert rc == 0

    schema_cols = list(ActivitiesSchema.columns)
    available = set(df.columns) | {
        "pipeline_version",
        "timestamp_utc",
        "lower_value",
        "upper_value",
        "activity_properties",
        "action_type",
        "properties_hash",
    }
    expected_head = [c for c in schema_cols if c in available]
    expected_tail = sorted(available - set(schema_cols))
    assert captured["col_order"] == expected_head + expected_tail


def test_run_chembl_streams_large_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure large datasets trigger chunked deterministic CSV writes."""

    chunk_size = cfg.io.csv_chunksize
    total_rows = chunk_size * 2 + 5
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n" + "\n".join(str(i) for i in range(total_rows)))

    monkeypatch.setattr(
        io,
        "read_ids",
        lambda *_, **__: (str(i) for i in range(total_rows)),
    )

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        data = list(ids)
        return pd.DataFrame(
            {
                "activity_id": data,
                "assay_chembl_id": [f"A{i}" for i in data],
                "molecule_chembl_id": [f"M{i}" for i in data],
                "standard_value": [1.0] * len(data),
                "standard_type": ["IC50"] * len(data),
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)
    monkeypatch.setattr(gad, "normalize_activities", lambda df: df)
    monkeypatch.setattr(gad, "add_pipeline_metadata", lambda df: df)
    monkeypatch.setattr(gad, "compute_activity_bounds", lambda df, cfg: df)
    monkeypatch.setattr(
        gad,
        "apply_activity_annotations",
        lambda df, action_cfg, properties_cfg: df,
    )

    class _Result:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    monkeypatch.setattr(
        gad,
        "validate_activities",
        lambda df, return_result: _Result(df),
    )
    monkeypatch.setattr(gad, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, object] = {}

    def fake_write_csv(
        df: Iterable[pd.DataFrame],
        output: Path,
        *,
        key_cols,
        col_order,
        chunksize,
        sort_chunksize,
        sep,
        encoding,
        cfg,
        **__,
    ) -> Path:
        frames = list(df)
        rows = sum(len(frame) for frame in frames)
        chunk_count = len(frames)
        captured.update(
            {
                "chunksize": chunksize,
                "sort_chunksize": sort_chunksize,
                "chunk_count": chunk_count,
                "rows": rows,
            }
        )
        return output

    monkeypatch.setattr(gad, "write_csv_chunks_deterministic", fake_write_csv)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")
    rc = gad.run_chembl(cfg, args)

    assert rc == 0
    assert captured["rows"] == total_rows
    assert captured["chunksize"] == chunk_size
    assert captured["sort_chunksize"] == chunk_size
    assert captured["chunk_count"] > 1


def test_run_chembl_workers_preserve_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("")

    cfg.activity.batch_size = 2
    cfg.activity.workers = 2

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(
        io, "read_ids", lambda *_, **__: (str(i) for i in range(1, 7))
    )

    delays = {("1", "2"): 0.05, ("3", "4"): 0.0, ("5", "6"): 0.0}

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        key = tuple(ids)
        time.sleep(delays.get(key, 0.0))
        return pd.DataFrame({"activity_id": list(ids)})

    monkeypatch.setattr(cl, "get_activities", fake_get)

    captured_chunks: list[list[str]] = []

    def fake_run_pipeline(*, fetcher, **kwargs):
        for chunk in fetcher():
            captured_chunks.append(list(chunk["activity_id"]))
        return 0

    monkeypatch.setattr(gad, "run_pipeline", fake_run_pipeline)

    rc = gad.run_chembl(cfg, args)

    assert rc == 0
    assert captured_chunks == [["1", "2"], ["3", "4"], ["5", "6"]]


def test_run_chembl_workers_chunk_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("")

    cfg.activity.batch_size = 2
    cfg.activity.workers = 2

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(
        io, "read_ids", lambda *_, **__: (str(i) for i in range(1, 7))
    )

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        items = list(ids)
        if items and items[0] == "3":
            raise ValueError("boom")
        return pd.DataFrame({"activity_id": items})

    monkeypatch.setattr(cl, "get_activities", fake_get)

    def fake_run_pipeline(*, fetcher, **kwargs):
        with pytest.raises(cli_utils.PipelineError):
            list(fetcher())
        return 1

    monkeypatch.setattr(gad, "run_pipeline", fake_run_pipeline)

    rc = gad.run_chembl(cfg, args)

    assert rc == 1


def test_run_chembl_workers_respect_rate_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("")

    cfg.activity.batch_size = 1
    cfg.activity.workers = 3
    cfg.sources.chembl.api.rps = 1
    cfg.sources.chembl.api.burst = 1

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(
        io, "read_ids", lambda *_, **__: (str(i) for i in range(3))
    )

    clock = SimpleNamespace(now=0.0, sleeps=[])

    def _monotonic() -> float:
        return clock.now

    def _sleep(delay: float) -> None:
        clock.sleeps.append(delay)
        clock.now += delay

    monkeypatch.setattr(rl, "time", SimpleNamespace(monotonic=_monotonic))
    monkeypatch.setattr(rl, "sleep", _sleep)

    with rl._limiters_lock:
        rl._limiters.clear()

    class DummyClient:
        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, *_args) -> bool:
            return False

        def request_json(self, url: str, *, cfg, timeout=None):
            limiter = rl.get_limiter("chembl", cfg.rps, cfg.burst)
            limiter.acquire()
            return {"activities": []}

    monkeypatch.setattr(gad, "ChemblClient", DummyClient)

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        client.request_json(f"chunk:{','.join(ids)}", cfg=cfg, timeout=timeout)
        return pd.DataFrame({"activity_id": list(ids)})

    monkeypatch.setattr(cl, "get_activities", fake_get)

    processed: list[list[str]] = []

    def fake_run_pipeline(*, fetcher, **kwargs):
        for chunk in fetcher():
            processed.append(list(chunk["activity_id"]))
        return 0

    monkeypatch.setattr(gad, "run_pipeline", fake_run_pipeline)

    rc = gad.run_chembl(cfg, args)

    assert rc == 0
    assert processed == [["0"], ["1"], ["2"]]
    assert clock.sleeps == pytest.approx([1.0, 1.0])


def test_run_chembl_integration_with_stubbed_request_json(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Integration test exercising ``run_chembl`` with the real pipeline."""

    input_csv = tmp_path / "activity.csv"
    output_csv = tmp_path / "output.csv"
    input_csv.write_text("activity_id\n1\n2\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv)

    monkeypatch.setattr(gad, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda _path: "stubbed")

    cfg.activity.batch_size = 2
    cfg.activity.workers = 1

    calls: list[tuple[str, ...]] = []

    def fake_request_json(self, url: str, *, cfg, timeout=None):
        parsed = urlparse(url)
        params = parse_qs(parsed.query)
        joined = params.get("activity_id__in", [""])[0]
        identifiers = tuple(filter(None, joined.split(",")))
        calls.append(identifiers)
        activities = [
            {
                "activity_id": identifier,
                "assay_chembl_id": f"ASSAY{identifier}",
                "molecule_chembl_id": f"CHEMBL{identifier}",
                "standard_type": "IC50",
                "standard_value": float(index + 1),
                "type": "IC50",
                "relation": "=",
                "value": None,
                "pchembl_value": None,
                "units": "nM",
            }
            for index, identifier in enumerate(identifiers)
        ]
        return {"activities": activities}

    monkeypatch.setattr(gad.ChemblClient, "request_json", fake_request_json)

    rc = gad.run_chembl(cfg, args)

    assert rc == 0
    assert calls == [("1", "2")]
    assert output_csv.exists()

    result = pd.read_csv(output_csv)
    assert list(result["activity_id"].astype(str)) == ["1", "2"]
    assert result["standard_value"].tolist() == [1.0, 2.0]
    assert result.empty is False
