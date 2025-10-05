from __future__ import annotations

import json
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
import yaml

import library.cli_utils as cli_utils
from library.config import Config
from library.pipelines.common import metadata as pipeline_metadata
from library.pipelines.tissue import chembl as tissue_chembl
from library.pipelines.tissue import pipeline as tissue_pipeline
from scripts import get_tissue_data


class _MemoryLogger:
    """In-memory logger collecting structured CLI events."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, Any]]] = []

    def debug(self, event: str, **fields: Any) -> None:
        self.events.append(("debug", event, dict(fields)))

    def info(self, event: str, **fields: Any) -> None:
        self.events.append(("info", event, dict(fields)))

    def warning(self, event: str, **fields: Any) -> None:
        self.events.append(("warning", event, dict(fields)))

    def error(self, event: str, **fields: Any) -> None:
        self.events.append(("error", event, dict(fields)))


def _load_responses(base: str, resource_dir: Path) -> dict[str, dict[str, Any]]:
    def _read(name: str) -> dict[str, Any]:
        return json.loads((resource_dir / name).read_text(encoding="utf-8"))

    return {
        f"{base}/tissue.json?format=json&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2": _read(
            "chembl_response_chunk0_page0.json"
        ),
        (
            f"{base}/tissue.json?format=json"
            "&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2&offset=1"
        ): _read("chembl_response_chunk0_page1.json"),
        f"{base}/tissue.json?format=json&tissue_chembl_id__in=CHEMBL999999&limit=1": _read(
            "chembl_response_chunk1_page0.json"
        ),
    }


class _StubChemblClient:
    def __init__(self, responses: dict[str, dict[str, Any]]) -> None:
        self._responses = responses
        self.calls: list[dict[str, Any]] = []
        self.closed = False

    def request_json(self, url: str, *, cfg: Any, timeout: float | None = None) -> dict[str, Any]:
        self.calls.append({"url": url, "timeout": timeout})
        if url not in self._responses:
            raise AssertionError(f"Unexpected request for {url}")
        return json.loads(json.dumps(self._responses[url]))

    def close(self) -> None:
        self.closed = True


@pytest.mark.e2e
def test_get_tissue_cli__smoke(tmp_path: Path, snapshot_resource: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch) -> None:
    resource_dir = snapshot_resource / "tissue"
    input_csv = tmp_path / "input_ids.csv"
    input_csv.write_text(
        (resource_dir / "input_ids.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    output_csv = tmp_path / "tissue_output.csv"

    cfg.io.output_dir = tmp_path
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.io.cache_dir.mkdir()
    cfg.io.exist_ok = True
    cfg.io.csv_encoding = "utf-8"
    cfg.tissue.batch_size = 2
    cfg.tissue.timeout = 15.0

    base = cfg.api.chembl_base.rstrip("/")
    responses = _load_responses(base, resource_dir)

    pipeline_metadata.get_timestamp_utc.cache_clear()
    pipeline_metadata.pipeline_metadata.cache_clear()
    monkeypatch.setattr(
        pipeline_metadata,
        "get_timestamp_utc",
        lambda: "2020-01-01T00:00:00+00:00",
    )
    pipeline_metadata.pipeline_metadata.cache_clear()

    frozen_now = datetime(2020, 1, 1, tzinfo=timezone.utc)

    def _frozen_datetime_now(tz: timezone | None = None) -> datetime:
        if tz is None:
            return frozen_now.replace(tzinfo=None)
        return frozen_now.astimezone(tz)

    monkeypatch.setattr(
        "library.io.metadata.datetime",
        SimpleNamespace(now=_frozen_datetime_now),
    )

    created_clients: list[_StubChemblClient] = []

    def _client_factory(*_args: Any, **_kwargs: Any) -> _StubChemblClient:
        client = _StubChemblClient(responses)
        created_clients.append(client)
        return client

    monkeypatch.setattr(tissue_pipeline, "ChemblClient", _client_factory)

    logger = _MemoryLogger()
    monkeypatch.setattr(get_tissue_data, "logger", logger)

    def _configure_logger(_log_cfg: Any) -> _MemoryLogger:
        return logger

    monkeypatch.setattr(get_tissue_data.cli, "configure_logger", _configure_logger)
    monkeypatch.setattr(get_tissue_data, "configure_logger", _configure_logger)
    monkeypatch.setattr(cli_utils.cli, "configure_logger", _configure_logger)

    @contextmanager
    def _fake_setup_cli_logging(script_name: str, log_cfg: Any, date_str: str | None = None, **_: Any):
        yield SimpleNamespace(log_cfg=log_cfg, console_stream=None)

    monkeypatch.setattr(get_tissue_data, "setup_cli_logging", _fake_setup_cli_logging)

    def _fake_apply_config_overrides(
        args: Any,
        parser: Any,
        config_path: Any,
        mapping: dict[str, str] | None = None,
        *,
        base_parser: Any | None = None,
    ) -> Config:
        args._config_metadata = SimpleNamespace(snapshot={"source": "test"})
        cfg.io.output_dir = tmp_path
        cfg.io.cache_dir = tmp_path / "cache"
        cfg.io.cache_dir.mkdir(exist_ok=True)
        cfg.io.exist_ok = True
        cfg.io.csv_encoding = "utf-8"
        cfg.tissue.batch_size = getattr(args, "batch_size", cfg.tissue.batch_size)
        cfg.tissue.limit = getattr(args, "limit", cfg.tissue.limit)
        cfg.tissue.offset = getattr(args, "offset", cfg.tissue.offset)
        if getattr(args, "timeout", None) is not None:
            cfg.tissue.timeout = float(args.timeout)
        return cfg

    monkeypatch.setattr(cli_utils, "apply_config_overrides", _fake_apply_config_overrides)

    argv = [
        "--config",
        str(tmp_path / "config.toml"),
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--batch-size",
        "2",
        "--mode",
        "chembl",
    ]

    before_first = len(logger.events)
    exit_code = get_tissue_data.main(argv)
    first_segment = logger.events[before_first:]

    assert exit_code == 0
    summary_events = [event for event in first_segment if event[1] == "tissue_pipeline_summary"]
    assert summary_events, "expected pipeline summary log"
    summary_payload = summary_events[-1][2]
    assert summary_payload["records"] == 3
    assert summary_payload["output"] == str(output_csv)
    assert float(summary_payload["duration"]) >= 0.0

    assert created_clients, "ChemblClient should be instantiated"
    for stub in created_clients:
        assert stub.closed is True
        assert [call["url"] for call in stub.calls] == list(responses.keys())

    output_snapshot = output_csv.read_text(encoding="utf-8")
    meta_path = output_csv.with_suffix(".csv.meta.yaml")
    assert meta_path.exists()
    meta_snapshot = meta_path.read_text(encoding="utf-8")
    meta = yaml.safe_load(meta_snapshot)
    assert meta["columns"] == tissue_chembl.TISSUE_COLUMN_ORDER

    before_second = len(logger.events)
    second_exit = get_tissue_data.main(argv)
    second_segment = logger.events[before_second:]

    assert second_exit == 0
    assert output_csv.read_text(encoding="utf-8") == output_snapshot
    assert meta_path.read_text(encoding="utf-8") == meta_snapshot

    second_summary = [event for event in second_segment if event[1] == "tissue_pipeline_summary"]
    assert second_summary
    second_payload = second_summary[-1][2]
    assert second_payload["records"] == 3
    assert second_payload["output"] == str(output_csv)
    assert float(second_payload["duration"]) >= 0.0

    exported = pd.read_csv(output_csv, dtype="string")
    assert list(exported.columns) == list(tissue_chembl.TISSUE_COLUMN_ORDER)
    assert list(exported["tissue_chembl_id"]) == [
        "CHEMBL2109249",
        "CHEMBL613507",
        "CHEMBL999999",
    ]
    for stub in created_clients:
        assert stub.closed is True
