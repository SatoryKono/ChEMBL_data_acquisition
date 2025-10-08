from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from library.config import Config
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
    input_csv.write_text("activity_id\nCHEMBL1\n", encoding="utf-8")
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
    input_csv.write_text("activity_id\nCHEMBL1\n", encoding="utf-8")
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
def test_main__config_limit_used_when_cli_omitted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    cfg.activity.limit = 4

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
    monkeypatch.setattr(
        get_activities.cli,
        "apply_config_overrides",
        lambda args, parser, config: cfg,
    )
    monkeypatch.setattr(get_activities, "ensure_dirs", lambda _cfg: None)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
        ]
    )

    assert exit_code == 0
    assert observed["limit"] == 4
    assert writer_calls == []
    assert (
        "info",
        "generated",
        {"count": 4, "output": str(output_csv)},
    ) in logger_stub.events
def test_parse_args__invalid_limit_error_message(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        get_activities.parse_args(["--limit", "-1"])

    captured = capsys.readouterr()
    assert "limit must be a non-negative integer" in captured.err


@pytest.mark.unit
def test_run__skip_existing_respects_flag(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "activities.csv"
    output_csv.write_text("activity_id\nOLD\n", encoding="utf-8")

    args = SimpleNamespace(
        limit=1,
        dry_run=False,
        skip_existing=True,
        force=False,
        output_csv=output_csv,
        input_csv=tmp_path / "input.csv",
    )

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)
    monkeypatch.setattr(
        get_activities,
        "get_activities",
        lambda *_args, **_kwargs: pytest.fail("get_activities should not be called"),
    )
    monkeypatch.setattr(
        get_activities,
        "_write_output",
        lambda *_args, **_kwargs: pytest.fail("_write_output should not be called"),
    )

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 0
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(output_csv)},
    ) in logger_stub.events
    assert output_csv.read_text(encoding="utf-8") == "activity_id\nOLD\n"


@pytest.mark.unit
def test_run__force_overrides_skip_existing(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "activities.csv"
    output_csv.write_text("activity_id\nOLD\n", encoding="utf-8")

    args = SimpleNamespace(
        limit=2,
        dry_run=False,
        skip_existing=True,
        force=True,
        output_csv=output_csv,
        input_csv=tmp_path / "input.csv",
    )

    monkeypatch.setattr(
        get_activities,
        "get_activities",
        lambda limit: [{"activity_id": idx + 1} for idx in range(limit)],
    )

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 0
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(output_csv)},
    ) not in logger_stub.events

    content = output_csv.read_text(encoding="utf-8")
    lines = content.splitlines()
    assert [lines[0].lstrip("\ufeff"), *lines[1:]] == ["activity_id", "1", "2"]

    meta_path = Path(f"{output_csv}.meta.yaml")
    assert meta_path.exists()

    leftover = list(output_csv.parent.glob("*.tmp"))
    leftover.extend(output_csv.parent.glob("*.tmp.meta.yaml"))
    assert not leftover
