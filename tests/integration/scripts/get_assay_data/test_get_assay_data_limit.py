from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library import io
from library.cli import utils as cli_utils
from library.config import Config
from scripts import get_assay_data as gas


class TrackingIterator:
    """Iterator that records how many identifiers were consumed."""

    def __init__(self, total: int) -> None:
        self._total = total
        self.consumed = 0

    def __iter__(self) -> TrackingIterator:
        return self

    def __next__(self) -> str:
        if self.consumed >= self._total:
            raise StopIteration
        self.consumed += 1
        return f"A{self.consumed}"


def test_run_chembl_limit_streams_ids_without_materialising(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """The limit option must avoid materialising identifiers into a list."""

    cfg.assay.limit = 2
    cfg.assay.batch_size = 1

    args = argparse.Namespace(
        input_csv=tmp_path / "assays.csv",
        output_csv=tmp_path / "out.csv",
    )
    args.input_csv.write_text("assay_chembl_id\nA1\n")

    iterator = TrackingIterator(total=10)
    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iterator)

    class DummyClient:
        def __enter__(self) -> DummyClient:
            return self

        def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> bool:
            return False

    monkeypatch.setattr(gas, "ChemblClient", lambda *a, **k: DummyClient())

    captured: dict[str, Any] = {}

    def fake_get_assays(
        ids_source, *, cfg, client, chunk_size, timeout
    ) -> pd.DataFrame:  # type: ignore[no-untyped-def]
        captured["ids_type"] = type(ids_source)
        ids_list = list(ids_source)
        captured["ids"] = ids_list
        return pd.DataFrame(
            {
                "assay_chembl_id": ids_list,
                "pipeline_version": ["1.0"] * len(ids_list),
                "timestamp_utc": ["now"] * len(ids_list),
            }
        )

    monkeypatch.setattr(gas.cl, "get_assays", fake_get_assays)
    monkeypatch.setattr(gas.ap, "postprocess_assays", lambda frame: frame)
    monkeypatch.setattr(gas, "normalize_assays", lambda frame: frame)
    monkeypatch.setattr(gas, "add_pipeline_metadata", lambda frame: frame)
    monkeypatch.setattr(gas, "analyze_table_quality", lambda *a, **k: None)

    class _Result:
        def __init__(self, data: pd.DataFrame) -> None:
            self.data = data
            self.failure_cases = pd.DataFrame()

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> _Result:
        return _Result(frame)

    monkeypatch.setattr(gas, "validate_assays", fake_validate)

    def fake_write(
        chunks,
        destination: Path,
        *,
        key_cols,
        col_order,
        chunksize,
        sort_chunksize,
        sep,
        encoding,
        cfg,
    ) -> Path:  # type: ignore[no-untyped-def]
        frames = list(chunks)
        assert frames
        values = frames[0]["assay_chembl_id"].tolist()
        destination.write_text("assay_chembl_id\n" + "\n".join(values))
        return destination

    monkeypatch.setattr(gas, "write_csv_chunks_deterministic", fake_write)
    monkeypatch.setattr(cli_utils, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(cli_utils, "file_sha256", lambda path: "deadbeef")

    log_calls: list[dict[str, Any]] = []

    def fake_info(
        event: str, *args: Any, extra: dict[str, Any] | None = None, **kwargs: Any
    ) -> None:
        log_calls.append({"event": event, "kwargs": kwargs, "extra": extra})

    monkeypatch.setattr(gas.logger, "info", fake_info)

    rc = gas.run_chembl(cfg, args)

    assert rc == 0
    assert iterator.consumed == cfg.assay.limit
    assert captured["ids_type"] is not list
    assert captured["ids"] == ["A1", "A2"]
    assert log_calls[-1]["event"] == "process_limit"
    assert log_calls[-1]["kwargs"] == {"limit": iterator.consumed}
