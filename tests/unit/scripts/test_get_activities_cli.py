from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

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

    def exception(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("exception", event, dict(data)))


def _make_cfg_stub(base: Path, *, limit: int | None = None, exist_ok: bool = True):
    io_cfg = SimpleNamespace(
        exist_ok=exist_ok,
        output_dir=base,
        cache_dir=base / "cache",
        csv_sep=",",
        csv_encoding="utf-8",
    )
    cfg_stub = SimpleNamespace(activity=SimpleNamespace(limit=limit), io=io_cfg)

    def _to_dict() -> dict[str, object]:
        return {
            "activity": {"limit": limit},
            "io": {
                "exist_ok": exist_ok,
                "output_dir": str(base),
                "cache_dir": str(base / "cache"),
                "csv_sep": ",",
                "csv_encoding": "utf-8",
            },
        }

    cfg_stub.to_dict = _to_dict  # type: ignore[attr-defined]
    return cfg_stub


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
    cfg_stub = _make_cfg_stub(tmp_path, limit=None)

    monkeypatch.setattr(
        get_activities.cli,
        "apply_config_overrides",
        lambda args, parser, config: cfg_stub,
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

    cfg_stub = _make_cfg_stub(tmp_path, limit=None)

    monkeypatch.setattr(
        get_activities.cli,
        "apply_config_overrides",
        lambda args, parser, config: cfg_stub,
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
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"

    cfg_stub = _make_cfg_stub(tmp_path, limit=4)

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
        lambda args, parser, config: cfg_stub,
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
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "activities.csv"
    output_csv.write_text("activity_id\nOLD\n", encoding="utf-8")

    cfg = _make_cfg_stub(tmp_path)

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
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_csv = tmp_path / "activities.csv"
    output_csv.write_text("activity_id\nOLD\n", encoding="utf-8")

    cfg = _make_cfg_stub(tmp_path)

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
def test_parse_args__only_expected_options_present() -> None:
    _parser, args, _log_cfg = get_activities.parse_args([])

    assert "chunk_size" not in vars(args)
    assert "column" not in vars(args)


@pytest.mark.unit
def test_run__default_limit_used_when_config_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_cfg_stub(tmp_path, limit=None)

    args = SimpleNamespace(
        limit=None,
        dry_run=False,
        skip_existing=False,
        force=False,
        output_csv=tmp_path / "activities.csv",
        input_csv=tmp_path / "input.csv",
    )

    observed: dict[str, int] = {}

    def _fake_get(limit: int) -> list[dict[str, int]]:
        observed["limit"] = limit
        return [{"activity_id": idx} for idx in range(limit)]

    monkeypatch.setattr(get_activities, "get_activities", _fake_get)
    monkeypatch.setattr(get_activities, "_write_output", lambda *_a, **_k: args.output_csv)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 0
    assert observed["limit"] == get_activities.DEFAULT_LIMIT


@pytest.mark.unit
def test_run__writes_output_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_cfg_stub(tmp_path, exist_ok=True)

    output_csv = tmp_path / "nested" / "activities.csv"
    args = SimpleNamespace(
        limit=2,
        dry_run=False,
        skip_existing=False,
        force=False,
        output_csv=output_csv,
        input_csv=tmp_path / "input.csv",
    )

    def _fake_get(limit: int) -> list[dict[str, int]]:
        return [{"activity_id": idx + 1} for idx in range(limit)]

    monkeypatch.setattr(get_activities, "get_activities", _fake_get)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 0
    assert output_csv.exists()
    content = output_csv.read_text(encoding="utf-8").splitlines()
    assert content == ["activity_id", "1", "2"]
    meta_path = Path(f"{output_csv}.meta.yaml")
    assert meta_path.exists()
    assert any(event[1] == "generated" for event in logger_stub.events)


@pytest.mark.unit
def test_run__fails_when_directory_missing_and_exist_ok_false(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _make_cfg_stub(tmp_path, exist_ok=False)

    output_dir = tmp_path / "missing" / "deeper"
    output_csv = output_dir / "activities.csv"

    args = SimpleNamespace(
        limit=1,
        dry_run=False,
        skip_existing=False,
        force=False,
        output_csv=output_csv,
        input_csv=tmp_path / "input.csv",
    )

    monkeypatch.setattr(
        get_activities,
        "get_activities",
        lambda *_a, **_k: pytest.fail("get_activities should not be called"),
    )

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    exit_code = get_activities.run(cfg, args)

    assert exit_code == 1
    assert not output_dir.exists()
    assert (
        "error",
        "output_directory_missing",
        {"directory": str(output_dir), "output": str(output_csv)},
    ) in logger_stub.events
