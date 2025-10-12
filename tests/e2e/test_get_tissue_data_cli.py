"""End-to-end tests for the ``get_tissue_data`` CLI."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from collections import deque
from collections.abc import Iterable
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

from library.common.run_context import RunContext, set_current
from library.pipelines.tissue import TISSUE_COLUMN_ORDER
from scripts import get_tissue_data


class MemoryLogger:
    """In-memory logger capturing structured events for assertions."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, Any]]] = []

    def bind(
        self, **_: Any
    ) -> MemoryLogger:  # pragma: no cover - interface compatibility
        return self

    def _store(
        self,
        level: str,
        event: str,
        *,
        extra: dict[str, Any] | None = None,
        **data: Any,
    ) -> None:
        payload = dict(extra or {})
        payload.update(data)
        self.events.append((level, event, payload))

    def log(self, level: str, event: str, **data: Any) -> None:
        extra = data.pop("extra", None)
        self._store(str(level).lower(), event, extra=extra, **data)

    def debug(
        self, event: str, **data: Any
    ) -> None:  # pragma: no cover - unused in assertions
        self._store("debug", event, extra=data.pop("extra", None), **data)

    def info(self, event: str, **data: Any) -> None:
        self._store("info", event, extra=data.pop("extra", None), **data)

    def warning(self, event: str, **data: Any) -> None:
        self._store("warning", event, extra=data.pop("extra", None), **data)

    warn = warning

    def error(self, event: str, **data: Any) -> None:
        self._store("error", event, extra=data.pop("extra", None), **data)

    def exception(
        self, event: str, **data: Any
    ) -> None:  # pragma: no cover - defensive fallback
        self._store("error", event, extra=data.pop("extra", None), **data)


class StubChemblClient:
    """Deterministic replacement for :class:`ChemblClient` used in tests."""

    def __init__(
        self, pages: Iterable[dict[str, Any]], calls: list[dict[str, Any]]
    ) -> None:
        self._pages = deque(copy.deepcopy(list(pages)))
        self._calls = calls

    def __enter__(self) -> StubChemblClient:  # pragma: no cover - interface parity
        return self

    def __exit__(
        self, exc_type, exc, tb
    ) -> bool:  # pragma: no cover - interface parity
        return False

    def close(self) -> None:  # pragma: no cover - interface parity
        return None

    @property
    def pages_remaining(self) -> int:
        return len(self._pages)

    def request_json(
        self, url: str, *, cfg: Any, timeout: float | None = None
    ) -> dict[str, Any]:
        del cfg  # the configuration is unused in the stub
        self._calls.append({"url": str(url), "timeout": timeout})
        if not self._pages:
            raise AssertionError(f"unexpected request: {url}")
        payload = self._pages.popleft()
        return copy.deepcopy(payload)


@pytest.mark.e2e
def test_get_tissue_data_cli__end_to_end(
    tmp_path: Path,
    cfg,
    snapshot_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resource_dir = snapshot_resource / "tissue_pipeline"
    page_files = [
        resource_dir / "chunk1_page1.json",
        resource_dir / "chunk1_page2.json",
        resource_dir / "chunk2_page1.json",
    ]
    base_pages = [json.loads(path.read_text(encoding="utf-8")) for path in page_files]

    recorded_calls: list[dict[str, Any]] = []
    created_clients: list[StubChemblClient] = []

    def _make_client(*_args: Any, **_kwargs: Any) -> StubChemblClient:
        client = StubChemblClient(base_pages, recorded_calls)
        created_clients.append(client)
        return client

    monkeypatch.setattr("library.pipelines.tissue.pipeline.ChemblClient", _make_client)

    logger = MemoryLogger()
    monkeypatch.setattr(get_tissue_data, "logger", logger)
    monkeypatch.setattr("library.common.log.logger", logger)
    monkeypatch.setattr("library.pipelines.tissue.pipeline.logger", logger)
    monkeypatch.setattr("library.pipelines.tissue.chembl.logger", logger)

    def fake_apply_config_overrides(
        args: Any,
        parser: Any,
        config_path: str | Path,
        mapping: dict[str, str] | None = None,
        *,
        base_parser: Any | None = None,
    ) -> Any:
        del parser, mapping, base_parser
        args._config_metadata = SimpleNamespace(snapshot={"config": str(config_path)})
        if hasattr(args, "input_csv"):
            args.input_csv = Path(args.input_csv)
            final_out = getattr(args, "final_out", None)
            if final_out is not None:
                final_path = Path(final_out)
                args.final_out = final_path
                args.output_csv = final_path
            else:
                legacy_out = getattr(args, "output_csv", None)
                if legacy_out is not None:
                    output_path = Path(legacy_out)
                    args.final_out = output_path
                    args.output_csv = output_path
        cfg.io.output_dir = tmp_path
        cfg.io.csv_sep = ","
        cfg.io.csv_encoding = "utf-8"
        cfg.tissue.batch_size = getattr(args, "batch_size", cfg.tissue.batch_size)
        cfg.tissue.column = getattr(args, "column", cfg.tissue.column)
        cfg.tissue.limit = getattr(args, "limit", cfg.tissue.limit)
        cfg.tissue.offset = getattr(args, "offset", cfg.tissue.offset)
        timeout = getattr(args, "timeout", None)
        if timeout is not None:
            cfg.tissue.timeout = timeout
        return cfg

    monkeypatch.setattr(
        "library.cli_utils.apply_config_overrides", fake_apply_config_overrides
    )
    monkeypatch.setattr("library.cli_utils.ensure_dirs", lambda _cfg: None)

    def fake_configure_logger(log_cfg: Any) -> MemoryLogger:
        logger.log(
            "debug", "configure_logger_called", log_level=getattr(log_cfg, "level", "")
        )
        set_current(
            RunContext(
                run_id=str(getattr(log_cfg, "run_id", "")),
                generated_at=str(getattr(log_cfg, "generated_at", "")),
            )
        )
        return logger

    monkeypatch.setattr("library.cli_utils.cli.configure_logger", fake_configure_logger)
    monkeypatch.setattr(get_tissue_data.cli, "configure_logger", fake_configure_logger)
    monkeypatch.setattr(get_tissue_data, "configure_logger", fake_configure_logger)

    @contextmanager
    def fake_setup_cli_logging(
        script_name: str, log_cfg: Any, date_str: str | None = None
    ):
        del script_name, date_str
        yield SimpleNamespace(log_cfg=log_cfg, console_stream=None)

    monkeypatch.setattr(get_tissue_data, "setup_cli_logging", fake_setup_cli_logging)

    artefacts: dict[str, object] = {}
    import library.pipelines.tissue.pipeline as tissue_pipeline

    original_write_csv = tissue_pipeline.write_csv_deterministic

    def _capture_write_csv(
        frame: pd.DataFrame,
        destination: Path,
        *,
        col_order: Iterable[str] | None = None,
        key_cols: Iterable[str] | None = None,
        chunksize: int | None = None,
        sort_chunksize: int | None = None,
        sep: str = ",",
        encoding: str = "utf-8",
        cfg=None,
        **kwargs,
    ) -> Path:
        written_path = original_write_csv(
            frame,
            destination,
            col_order=col_order,
            key_cols=key_cols,
            chunksize=chunksize,
            sort_chunksize=sort_chunksize,
            sep=sep,
            encoding=encoding,
            cfg=cfg,
            **kwargs,
        )
        target = Path(destination)
        if target == output_csv:
            dataset = frame.copy()
            quality_df = pd.DataFrame(
                [
                    {"metric": "rows_total", "value": int(dataset.shape[0])}
                ]
            )
            correlation_df = pd.DataFrame(
                [
                    {
                        "column_a": "tissue_chembl_id",
                        "column_b": "pref_name",
                        "correlation": 0.0,
                    }
                ]
            )
            quality_path = target.with_name(
                f"{target.stem}_quality_report_table.csv"
            )
            correlation_path = target.with_name(
                f"{target.stem}_data_correlation_report_table.csv"
            )
            quality_df.to_csv(quality_path, index=False, sep=sep, encoding=encoding)
            correlation_df.to_csv(
                correlation_path, index=False, sep=sep, encoding=encoding
            )
            artefacts.update(
                {
                    "dataset": dataset,
                    "quality": quality_df,
                    "correlation": correlation_df,
                    "paths": {
                        "dataset": written_path,
                        "quality": quality_path,
                        "correlation": correlation_path,
                    },
                }
            )
        return written_path

    monkeypatch.setattr(
        "library.pipelines.tissue.pipeline.write_csv_deterministic",
        _capture_write_csv,
    )

    # Ensure metadata timestamps are deterministic across repeated invocations.
    monkeypatch.setenv("CHEMBL_DA_RUN_ID", "tissue-e2e-test")

    input_csv = tmp_path / "inputs" / "tissues.csv"
    input_csv.parent.mkdir(parents=True, exist_ok=True)
    input_csv.write_text(
        "tissue_chembl_id\nCHEMBLT1\nCHEMBLT2\nCHEMBLT3\n",
        encoding="utf-8",
    )

    output_csv = tmp_path / "outputs" / "tissue.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    config_path = tmp_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n", encoding="utf-8")

    def _invoke(argv: list[str]) -> int:
        args = [str(part) for part in argv]
        monkeypatch.setattr(sys, "argv", ["get_tissue_data", *args])
        return get_tissue_data.main(args)

    base_argv = [
        "--config",
        str(config_path),
        "--input",
        str(input_csv),
        "--final-out",
        str(output_csv),
        "--batch-size",
        "2",
    ]

    exit_code_first = _invoke(base_argv)
    assert exit_code_first == 0

    assert created_clients, "expected ChemblClient to be instantiated"
    assert recorded_calls, "expected stub client to be exercised"
    first_run_calls = recorded_calls[: len(base_pages)]
    assert first_run_calls[0]["url"].endswith("CHEMBLT1,CHEMBLT2&limit=2")
    assert "CHEMBLT3" in first_run_calls[-1]["url"]

    fetch_event = next(
        (
            payload
            for level, event, payload in logger.events
            if event == "tissue_fetch_start"
        ),
        None,
    )
    assert fetch_event is not None
    assert fetch_event["requested"] == 3

    start_event = next(
        (event for event in logger.events if event[1] == "tissue_pipeline_start"), None
    )
    assert start_event is not None
    assert start_event[2]["input"] == str(input_csv)
    assert start_event[2]["output"] == str(output_csv)

    summary_event = next(
        (event for event in logger.events if event[1] == "tissue_pipeline_summary"),
        None,
    )
    assert summary_event is not None
    assert summary_event[2]["records"] == 3
    assert summary_event[2]["duration"] >= 0

    warning_events = {event for level, event, _ in logger.events if level == "warning"}
    assert "tissue_missing_identifiers" in warning_events
    assert "tissue_missing_identifiers_summary" in warning_events

    output_df = pd.read_csv(output_csv, sep=cfg.io.csv_sep, dtype="string")
    assert list(output_df.columns) == TISSUE_COLUMN_ORDER
    assert output_df["tissue_chembl_id"].tolist() == [
        "CHEMBLT1",
        "CHEMBLT2",
        "CHEMBLT3",
    ]
    assert (
        output_df.loc[output_df["tissue_chembl_id"] == "CHEMBLT3", "pref_name"]
        .isna()
        .all()
    )

    failure_path = output_csv.with_name(f"{output_csv.stem}_validation_failures.csv")
    assert not failure_path.exists()
    recorded_paths = artefacts.get("paths")
    assert recorded_paths is not None
    quality_path = recorded_paths["quality"]
    correlation_path = recorded_paths["correlation"]
    assert Path(recorded_paths["dataset"]) == output_csv
    assert quality_path.exists()
    assert correlation_path.exists()

    produced_files = {path.name for path in output_csv.parent.glob("*.csv")}
    assert produced_files == {
        output_csv.name,
        quality_path.name,
        correlation_path.name,
    }
    assert not Path(f"{failure_path}.meta.yaml").exists()
    assert not list(output_csv.parent.glob("*.meta.yaml"))

    csv_hash_first = hashlib.sha256(output_csv.read_bytes()).hexdigest()
    quality_hash_first = hashlib.sha256(quality_path.read_bytes()).hexdigest()
    correlation_hash_first = hashlib.sha256(correlation_path.read_bytes()).hexdigest()
    first_event_count = len(logger.events)

    exit_code_second = _invoke(base_argv)
    assert exit_code_second == 0

    assert len(created_clients) == 2
    assert all(client.pages_remaining == 0 for client in created_clients)
    assert len(recorded_calls) == len(base_pages) * 2

    csv_hash_second = hashlib.sha256(output_csv.read_bytes()).hexdigest()
    assert csv_hash_second == csv_hash_first

    updated_paths = artefacts.get("paths")
    assert updated_paths is not None
    assert Path(updated_paths["dataset"]) == output_csv
    assert Path(updated_paths["quality"]) == quality_path
    assert Path(updated_paths["correlation"]) == correlation_path
    assert hashlib.sha256(output_csv.read_bytes()).hexdigest() == csv_hash_first
    assert hashlib.sha256(quality_path.read_bytes()).hexdigest() == quality_hash_first
    assert hashlib.sha256(correlation_path.read_bytes()).hexdigest() == correlation_hash_first

    assert not failure_path.exists()
    assert not Path(f"{failure_path}.meta.yaml").exists()

    new_events = logger.events[first_event_count:]
    assert any(event == "tissue_pipeline_start" for _, event, _ in new_events)
    assert any(
        event == "tissue_pipeline_summary"
        and data.get("records") == 3
        and data.get("duration", 0) >= 0
        for _, event, data in new_events
    )
