from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
import pandas as pd

from library.utils.cli_tools import get_activities
from library.config import Config


class _LoggerStub:
    """Capture structured logging events emitted by the CLI."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def bind(self, **_: object) -> "_LoggerStub":  # pragma: no cover - interface parity
        return self

    def info(self, event: str, **data: object) -> None:
        self.events.append(("info", event, dict(data)))

    def warning(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("warning", event, dict(data)))

    def error(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("error", event, dict(data)))


@pytest.mark.unit
def test_main__limit_forwarded_to_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    observed: dict[str, int] = {}

    def _fake_get_activities(limit: int) -> list[dict[str, int]]:
        observed["limit"] = limit
        return [{"activity_id": idx} for idx in range(limit)]

    monkeypatch.setattr(
        "library.pipelines.activity.get_activities",
        _fake_get_activities,
    )
    monkeypatch.setattr(get_activities, "get_activities", _fake_get_activities)

    writer_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []
    original_prepare = get_activities.cli.prepare_io_paths

    def _prepare_io_paths(args: object, *, output_stem: str | None = None) -> None:
        original_prepare(args, output_stem=output_stem)
        setattr(args, "writer", lambda *a, **kw: writer_calls.append((a, kw)))

    monkeypatch.setattr(get_activities.cli, "prepare_io_paths", _prepare_io_paths)
    monkeypatch.setattr(get_activities.cli, "configure_logger", lambda *_a, **_k: None)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "7",
        ]
    )

    assert exit_code == 0
    assert observed["limit"] == 7
    assert writer_calls == []
    assert (
        "info",
        "generated",
        {"count": 7, "output": str(output_csv)},
    ) in logger_stub.events


@pytest.mark.unit
def test_main__dry_run_skips_fetch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    called: dict[str, bool] = {"value": False}

    def _unexpected_call(limit: int) -> list[dict[str, int]]:  # pragma: no cover - defensive
        called["value"] = True
        return [{"activity_id": limit}]

    monkeypatch.setattr(
        "library.pipelines.activity.get_activities",
        _unexpected_call,
    )
    monkeypatch.setattr(get_activities, "get_activities", _unexpected_call)

    writer_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []
    original_prepare = get_activities.cli.prepare_io_paths

    def _prepare_io_paths(args: object, *, output_stem: str | None = None) -> None:
        original_prepare(args, output_stem=output_stem)
        setattr(args, "writer", lambda *a, **kw: writer_calls.append((a, kw)))

    monkeypatch.setattr(get_activities.cli, "prepare_io_paths", _prepare_io_paths)
    monkeypatch.setattr(get_activities.cli, "configure_logger", lambda *_a, **_k: None)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "5",
            "--dry-run",
        ]
    )

    assert exit_code == 0
    assert called["value"] is False
    assert writer_calls == []
    assert ("info", "dry_run", {"limit": 5}) in logger_stub.events


@pytest.mark.unit
def test_run__validation_error_returns_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = Config()
    output_csv = tmp_path / "activities.csv"
    args = SimpleNamespace(
        limit=3,
        dry_run=False,
        output_csv=output_csv,
        input_csv=None,
    )

    monkeypatch.setattr(
        get_activities,
        "get_activities",
        lambda limit: [{"activity_id": idx} for idx in range(limit)],
    )

    def _fail_validation(_: pd.DataFrame, *, return_result: bool) -> object:
        raise ValueError("broken")

    monkeypatch.setattr(get_activities, "validate_activities", _fail_validation)

    table_quality_calls: list[object] = []

    def _noop_quality(*_args: object, **_kwargs: object):
        table_quality_calls.append(None)
        return lambda *_a, **_k: None

    monkeypatch.setattr(get_activities, "build_table_quality_hook", _noop_quality)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 1
    assert table_quality_calls == []


@pytest.mark.unit
def test_run__table_quality_invoked_after_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = Config()
    output_csv = tmp_path / "activities.csv"
    args = SimpleNamespace(
        limit=2,
        dry_run=False,
        output_csv=output_csv,
        input_csv=None,
    )

    records = [
        {
            "activity_id": 0,
            "molecule_chembl_id": "CHEMBL0",
            "assay_chembl_id": "ASSAY0",
            "standard_value": 1.0,
        },
        {
            "activity_id": 1,
            "molecule_chembl_id": "CHEMBL1",
            "assay_chembl_id": "ASSAY1",
            "standard_value": 2.0,
        },
    ]

    monkeypatch.setattr(get_activities, "get_activities", lambda limit: records[:limit])

    observed: dict[str, object] = {}

    def _build_quality(cfg_section, *, table_name, destination):
        observed["cfg"] = cfg_section
        observed["table_name"] = table_name
        observed["destination"] = destination

        def _hook(subject: object) -> None:
            observed["subject"] = subject

        return _hook

    monkeypatch.setattr(get_activities, "build_table_quality_hook", _build_quality)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 0
    assert output_csv.exists()
    assert observed["cfg"] == cfg.system.doc_quality
    assert observed["subject"] == output_csv
    assert observed["table_name"] == output_csv.with_suffix("")
    assert observed["destination"] == output_csv.parent
