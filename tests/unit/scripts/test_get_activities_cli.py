from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from library.utils.cli_tools import get_activities


class _LoggerStub:
    """Capture structured logging events emitted by the CLI."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def bind(self, **_: object) -> "_LoggerStub":  # pragma: no cover - interface parity
        return self

    def info(self, event: str, **data: object) -> None:
        self.events.append(("info", event, dict(data)))

    def warning(
        self, event: str, **data: object
    ) -> None:  # pragma: no cover - defensive
        self.events.append(("warning", event, dict(data)))

    def error(self, event: str, **data: object) -> None:  # pragma: no cover - defensive
        self.events.append(("error", event, dict(data)))

    def exception(
        self, event: str, **data: object
    ) -> None:  # pragma: no cover - defensive
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

    def _unexpected_call(
        limit: int,
    ) -> list[dict[str, int]]:  # pragma: no cover - defensive
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
def test_main__skip_existing_short_circuit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    output_csv.write_text("existing\n", encoding="utf-8")

    def _fake_get_activities(limit: int) -> list[dict[str, int]]:
        return [{"activity_id": value} for value in range(limit)]

    monkeypatch.setattr(
        "library.pipelines.activity.get_activities",
        _fake_get_activities,
    )
    monkeypatch.setattr(get_activities, "get_activities", _fake_get_activities)

    original_prepare = get_activities.cli.prepare_io_paths

    def _prepare_io_paths(args: object, *, output_stem: str | None = None) -> None:
        original_prepare(args, output_stem=output_stem)

    monkeypatch.setattr(get_activities.cli, "prepare_io_paths", _prepare_io_paths)
    monkeypatch.setattr(get_activities.cli, "configure_logger", lambda *_a, **_k: None)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    called = {"value": False}

    cfg_stub = _make_cfg_stub(tmp_path)
    monkeypatch.setattr(
        get_activities.cli,
        "apply_config_overrides",
        lambda args, parser, config: cfg_stub,
    )
    monkeypatch.setattr(get_activities, "ensure_dirs", lambda _cfg: None)

    def _unexpected_write(
        *_: object, **__: object
    ) -> Path:  # pragma: no cover - defensive
        called["value"] = True
        raise AssertionError("_write_output should not be invoked when skipping")

    monkeypatch.setattr(get_activities, "_write_output", _unexpected_write)

    exit_code = get_activities.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "3",
            "--skip-existing",
        ]
    )

    assert exit_code == 0
    assert called["value"] is False
    assert output_csv.read_text(encoding="utf-8") == "existing\n"
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(output_csv)},
    ) in logger_stub.events


@pytest.mark.unit
def test_main__skip_existing_with_force_writes_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    output_csv.write_text("existing\n", encoding="utf-8")

    def _fake_get_activities(limit: int) -> list[dict[str, int]]:
        return [{"activity_id": index} for index in range(limit)]

    monkeypatch.setattr(
        "library.pipelines.activity.get_activities",
        _fake_get_activities,
    )
    monkeypatch.setattr(get_activities, "get_activities", _fake_get_activities)

    original_prepare = get_activities.cli.prepare_io_paths

    def _prepare_io_paths(args: object, *, output_stem: str | None = None) -> None:
        original_prepare(args, output_stem=output_stem)

    monkeypatch.setattr(get_activities.cli, "prepare_io_paths", _prepare_io_paths)
    monkeypatch.setattr(get_activities.cli, "configure_logger", lambda *_a, **_k: None)

    logger_stub = _LoggerStub()
    monkeypatch.setattr(get_activities, "logger", logger_stub)

    cfg_stub = _make_cfg_stub(tmp_path, limit=None)
    monkeypatch.setattr(
        get_activities.cli,
        "apply_config_overrides",
        lambda args, parser, config: cfg_stub,
    )
    monkeypatch.setattr(get_activities, "ensure_dirs", lambda _cfg: None)

    write_calls: list[Path] = []

    def _write_csv(frame, path, *, cfg):  # type: ignore[override]
        path_obj = Path(path)
        write_calls.append(path_obj)
        frame.to_csv(path_obj, index=False)
        meta_path = path_obj.with_suffix(path_obj.suffix + ".meta.yaml")
        meta_path.write_text(
            yaml.safe_dump(
                {
                    "columns": list(frame.columns),
                    "dtypes": {
                        column: str(dtype) for column, dtype in frame.dtypes.items()
                    },
                },
                sort_keys=False,
            ),
            encoding="utf-8",
        )
        return path_obj

    monkeypatch.setattr(get_activities.io, "write_csv", _write_csv)

    exit_code = get_activities.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--limit",
            "2",
            "--skip-existing",
            "--force",
        ]
    )

    assert exit_code == 0
    assert write_calls, "write_csv should be invoked"
    assert write_calls[0] != output_csv
    assert "activity_id" in output_csv.read_text(encoding="utf-8")
    assert any(event[1] == "activity_generated" for event in logger_stub.events)


@pytest.mark.unit
def test_write_output__creates_meta_sidecar_and_removes_temporaries(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "activities.csv"
    frame = pd.DataFrame(
        {
            "activity_id": pd.Series([1, 2], dtype="Int64"),
            "name": pd.Series(["alpha", "beta"], dtype="string"),
        }
    )

    cfg_stub = _make_cfg_stub(tmp_path)

    result = get_activities._write_output(frame, output_path, cfg=cfg_stub)

    assert result == output_path
    assert output_path.exists()

    meta_path = output_path.with_suffix(output_path.suffix + ".meta.yaml")
    assert meta_path.exists()

    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8")) or {}

    assert meta.get("columns") == list(frame.columns)
    assert meta.get("dtypes") == {
        column: str(dtype) for column, dtype in frame.dtypes.items()
    }

    leftover_temporaries = [
        child.name
        for child in tmp_path.iterdir()
        if child.is_file()
        and (
            child.name.startswith(f".{output_path.name}.")
            or child.name.endswith(".tmp")
        )
    ]
    assert leftover_temporaries == []


@pytest.mark.unit
def test_write_output__fails_on_metadata_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output_path = tmp_path / "activities.csv"
    frame = pd.DataFrame(
        {
            "activity_id": pd.Series([1, 2], dtype="Int64"),
            "name": pd.Series(["alpha", "beta"], dtype="string"),
        }
    )

    cfg_stub = _make_cfg_stub(tmp_path)

    def _write_csv(frame: pd.DataFrame, path: Path, *, cfg):  # type: ignore[override]
        path_obj = Path(path)
        frame.to_csv(path_obj, index=False)
        meta_path = path_obj.with_suffix(path_obj.suffix + ".meta.yaml")
        meta_path.write_text(
            yaml.safe_dump(
                {"columns": ["wrong"], "dtypes": {"wrong": "string"}},
                sort_keys=False,
            ),
            encoding="utf-8",
        )
        return path_obj

    monkeypatch.setattr(get_activities.io, "write_csv", _write_csv)

    with pytest.raises(ValueError, match="metadata sidecar schema mismatch"):
        get_activities._write_output(frame.copy(), output_path, cfg=cfg_stub)

    assert list(tmp_path.glob("*.csv")) == []
    assert list(tmp_path.glob("*.meta.yaml")) == []
