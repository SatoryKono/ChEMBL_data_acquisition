"""CLI tests for :mod:`library.utils.cli_tools.pipeline_targets_main`."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library.config import Config
from library.pipeline_targets import PipelineResult
from library.utils.cli_tools import pipeline_targets_main as cli


class _DummyLogger:
    def __init__(
        self,
        context: dict[str, Any] | None = None,
        storage: list[tuple[str, dict[str, Any]]] | None = None,
    ) -> None:
        self._context: dict[str, Any] = context or {}
        self._records: list[tuple[str, dict[str, Any]]] = (
            storage if storage is not None else []
        )

    def bind(self, **ctx: Any) -> _DummyLogger:  # pragma: no cover - trivial
        merged = {**self._context, **ctx}
        return _DummyLogger(merged, self._records)

    def info(
        self, event: str, *args: Any, **kwargs: Any
    ) -> None:  # pragma: no cover - trivial
        record = {**self._context, **kwargs}
        self._records.append((event, record))

    @property
    def records(self) -> list[tuple[str, dict[str, Any]]]:  # pragma: no cover - trivial
        return list(self._records)


def test_cli_forwards_batch_size(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--batch-size",
        "25",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 25
    assert captured["chunks"] == [["CHEMBL1"]]
    assert captured["written_path"] == output_csv
    assert list(captured["written_df"]["target_chembl_id"]) == ["CHEMBL1"]
    assert any(
        event == "pipeline_start" and rec.get("stage") == "pipeline"
        for event, rec in dummy_logger.records
    )
    assert any(
        event == "pipeline_done"
        and rec.get("stage") == "pipeline"
        and rec.get("exit_code") == 0
        for event, rec in dummy_logger.records
    )


def test_cli_limit_restricts_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\nCHEMBL3\n", encoding="utf8"
    )
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--limit",
        "2",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 100
    assert captured["chunks"] == [["CHEMBL1", "CHEMBL2"]]
    assert captured["written_path"] == output_csv


def test_cli_limit_allows_zero(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\n", encoding="utf8"
    )
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(
            chembl=pd.DataFrame({"target_chembl_id": pd.Series(dtype="string")})
        )

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--limit",
        "0",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 100
    assert captured["chunks"] == []
    assert captured["written_path"] == output_csv
    assert captured["written_df"].empty


def test_cli_does_not_print_config_when_flag_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        next(chunk_iterator())
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        if not isinstance(data, pd.DataFrame):
            list(data)
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)

    def fail_print_config(_: Config) -> None:
        raise AssertionError("print_config should not be called")

    monkeypatch.setattr(cli, "print_config", fail_print_config)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out == ""


def test_cached_chembl_fetch_uses_chunk_concatenation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    chunks = iter(
        [
            ["CHEMBL1", "CHEMBL2"],
            ["CHEMBL3"],
        ]
    )

    recorded: dict[str, Any] = {}
    original_concat = pd.concat

    def spy_concat(objs: Iterable[pd.DataFrame], **kwargs: Any) -> pd.DataFrame:
        recorded["type"] = type(objs)
        frames = list(objs)
        recorded["sizes"] = [len(frame) for frame in frames]
        return original_concat(frames, **kwargs)

    monkeypatch.setattr(pd, "concat", spy_concat)

    df = cli._cached_chembl_fetch(chunks, Config())

    assert list(df["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    assert list(df["source"]) == ["chembl", "chembl", "chembl"]
    assert recorded["type"].__name__ == "chain"
    assert recorded["sizes"] == [2, 1]


def test_write_outputs_streams_large_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = Config()
    options = cli.PipelineConfig(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "out.csv",
        chunk_size=50,
        batch_size=50,
    )

    def frame_iterator() -> Iterable[pd.DataFrame]:
        for idx in range(0, 500, 75):
            upper = min(idx + 75, 500)
            ids = [f"CHEMBL{value}" for value in range(idx, upper)]
            yield pd.DataFrame({"target_chembl_id": ids})

    recorded: dict[str, Any] = {}

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Path:
        recorded["path"] = Path(path)
        recorded["is_generator"] = not isinstance(data, pd.DataFrame)
        frames = (
            [chunk.copy() for chunk in data]
            if recorded["is_generator"]
            else [data.copy()]
        )
        recorded["chunk_sizes"] = [len(frame) for frame in frames]
        return Path(path)

    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    output = cli._write_outputs(cfg, options, frame_iterator())

    assert output == options.output_csv
    assert recorded["path"] == options.output_csv
    assert recorded["is_generator"] is True
    assert recorded["chunk_sizes"] == [75, 75, 75, 75, 75, 75, 50]
