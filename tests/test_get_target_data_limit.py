"""Limit handling tests for :mod:`scripts.get_target_data`."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.config import Config
from scripts import get_target_data as gtd


class DummyClient:
    """Minimal context manager replacing :class:`~library.clients.ChemblClient`."""

    def __init__(
        self, *args: object, **kwargs: object
    ) -> None:  # pragma: no cover - trivial
        pass

    def __enter__(self) -> object:  # pragma: no cover - trivial
        return object()

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial
        return False


class DummyLogger:
    """Capture structured logging calls for assertions."""

    def __init__(self) -> None:
        self.calls: list[tuple[str, str, dict[str, Any]]] = []

    def info(self, event: str, **kwargs: Any) -> None:
        self.calls.append(("info", event, kwargs))

    def warning(self, event: str, **kwargs: Any) -> None:  # pragma: no cover - unused
        self.calls.append(("warning", event, kwargs))

    def debug(self, event: str, **kwargs: Any) -> None:  # pragma: no cover - unused
        self.calls.append(("debug", event, kwargs))

    def error(self, event: str, **kwargs: Any) -> None:  # pragma: no cover - unexpected
        raise AssertionError(f"Unexpected error log: {event!r} {kwargs!r}")


def test_main_limit_zero_skips_pipeline(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """CLI should short-circuit when ``--limit 0`` is provided."""

    recorded: list[tuple[str, dict[str, object]]] = []

    def capture_info(event: str, **kwargs: object) -> None:
        recorded.append((event, kwargs))

    monkeypatch.setattr(gtd.logger, "info", capture_info)

    def fail_run_cli_command(*_: object, **__: object) -> int:
        pytest.fail("run_cli_command must not execute when limit is zero")

    monkeypatch.setattr(gtd, "run_cli_command", fail_run_cli_command)

    output_csv = tmp_path / "targets.csv"

    exit_code = gtd.main(["chembl", "--limit", "0", "--output", str(output_csv)])

    assert exit_code == 0
    assert recorded == [("pipeline_skip_limit", {"limit": 0})]
    assert not output_csv.exists()
    assert not Path(f"{output_csv}.meta.yaml").exists()


def test_run_chembl_limit_uses_generator(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    """``run_chembl`` limits identifiers lazily and logs the processed count."""

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\nCHEMBL3\n",
        encoding=cfg.io.csv_encoding,
    )
    output_csv = tmp_path / "out.csv"

    config = cfg
    config.target.chembl.column = "target_chembl_id"
    config.target.chembl.limit = 2

    source_ids = ("CHEMBL1", "CHEMBL2", "CHEMBL3")

    def fake_read_ids(
        path: Path,
        *,
        column: str,
        cfg: object,
        **_: object,
    ) -> Any:
        assert path == input_csv
        assert column == config.target.chembl.column

        def generator() -> Any:
            for value in source_ids:
                yield value

        return generator()

    consumed: list[str] = []

    def fake_get_targets(ids_iter: Any, **__: object) -> pd.DataFrame:
        if isinstance(ids_iter, list):
            assert len(ids_iter) <= config.target.chembl.chunk_size
        for value in ids_iter:
            consumed.append(value)
        return pd.DataFrame({"target_chembl_id": consumed})

    dummy_logger = DummyLogger()
    monkeypatch.setattr(gtd, "logger", dummy_logger)
    monkeypatch.setattr(gtd.io, "read_ids", fake_read_ids)
    monkeypatch.setattr(gtd.cl, "get_targets", fake_get_targets)
    monkeypatch.setattr(gtd, "ChemblClient", DummyClient)
    monkeypatch.setattr(gtd, "normalize_targets", lambda df: df)
    monkeypatch.setattr(gtd, "add_pipeline_metadata", lambda df: df)
    monkeypatch.setattr(
        gtd, "_prepare_targets_for_schema", lambda df: (df, set(), set())
    )
    monkeypatch.setattr(
        gtd.TargetsSchema, "validate", staticmethod(lambda df, lazy=True: df)
    )
    def fake_write_csv(df: pd.DataFrame, path: Path, **__: object) -> Path:
        path.write_text("target_chembl_id\n", encoding="utf-8")
        return path

    monkeypatch.setattr(gtd.io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)

    args = argparse.Namespace(input_csv=input_csv, output_csv=output_csv, offset=0)

    exit_code = gtd.run_chembl(config, args)

    assert exit_code == 0
    assert consumed == ["CHEMBL1", "CHEMBL2"]
    assert any(
        level == "info"
        and event == "process_limit"
        and call_kwargs.get("limit") == len(consumed)
        for level, event, call_kwargs in dummy_logger.calls
    )
