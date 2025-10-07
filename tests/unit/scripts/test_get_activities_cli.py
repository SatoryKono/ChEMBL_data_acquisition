from __future__ import annotations

from pathlib import Path

import pytest

from library.utils.cli_tools import get_activities


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
