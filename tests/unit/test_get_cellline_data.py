import argparse
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.config import Config
from library.pipelines.cellline import CellLinePipelineOptions, CellLinePipelineResult
from scripts import get_cellline_data


class _MemoryLogger:
    """Capture structured log events emitted by the cell line CLI."""

    def __init__(self) -> None:
        self.events: list[tuple[str, str, dict[str, object]]] = []

    def info(self, event: str, **payload: object) -> None:
        self.events.append(("info", event, dict(payload)))

    def warning(self, event: str, **payload: object) -> None:
        self.events.append(("warning", event, dict(payload)))

    def error(self, event: str, **payload: object) -> None:
        self.events.append(("error", event, dict(payload)))


@pytest.fixture()
def logger_stub(monkeypatch: pytest.MonkeyPatch) -> _MemoryLogger:
    logger = _MemoryLogger()
    monkeypatch.setattr(get_cellline_data, "logger", logger)
    return logger


@pytest.fixture()
def minimal_args(tmp_path: Path) -> argparse.Namespace:
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("cell_chembl_id\nCHEMBL1\n", encoding="utf-8")
    output_csv = tmp_path / "output.csv"
    return argparse.Namespace(
        input_csv=input_csv,
        output_csv=output_csv,
        skip_existing=False,
        force=False,
        offset=0,
        limit=None,
        mode="chembl",
        batch_size=5,
        timeout=None,
        column="cell_chembl_id",
    )


def test_run__unsupported_mode(
    cfg: Config, minimal_args: argparse.Namespace, logger_stub: _MemoryLogger
) -> None:
    minimal_args.mode = "invalid"

    exit_code = get_cellline_data.run(cfg, minimal_args)

    assert exit_code == 1
    assert (
        "error",
        "cellline_unsupported_mode",
        {"mode": "invalid"},
    ) in logger_stub.events


def test_run__skip_existing_returns_zero(
    cfg: Config,
    minimal_args: argparse.Namespace,
    logger_stub: _MemoryLogger,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    minimal_args.skip_existing = True
    minimal_args.force = False
    minimal_args.output_csv.write_text("existing", encoding="utf-8")

    called = False

    def fake_run_pipeline(*_: object, **__: object) -> CellLinePipelineResult:
        nonlocal called
        called = True
        return CellLinePipelineResult(
            exit_code=0,
            records=0,
            duration=0.0,
            output_path=minimal_args.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=False,
        )

    monkeypatch.setattr(get_cellline_data, "run_cellline_pipeline", fake_run_pipeline)

    exit_code = get_cellline_data.run(cfg, minimal_args)

    assert exit_code == 0
    assert not called
    assert (
        "info",
        "pipeline_skip_existing",
        {"output": str(minimal_args.output_csv)},
    ) in logger_stub.events


def test_run__invokes_pipeline_with_expected_options(
    cfg: Config,
    minimal_args: argparse.Namespace,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_run_pipeline(
        cfg_arg: Config, options: CellLinePipelineOptions
    ) -> CellLinePipelineResult:
        captured["cfg"] = cfg_arg
        captured["options"] = options
        return CellLinePipelineResult(
            exit_code=0,
            records=2,
            duration=1.0,
            output_path=options.output_csv,
            failure_path=None,
            failures=0,
            missing_ids=(),
            written=True,
        )

    monkeypatch.setattr(get_cellline_data, "run_cellline_pipeline", fake_run_pipeline)

    minimal_args.limit = 5
    minimal_args.offset = 1
    minimal_args.batch_size = 7
    minimal_args.timeout = 12.5

    exit_code = get_cellline_data.run(cfg, minimal_args)

    assert exit_code == 0
    assert captured["cfg"] is cfg
    options = captured["options"]
    assert isinstance(options, CellLinePipelineOptions)
    assert options.limit == 5
    assert options.offset == 1
    assert options.batch_size == 7
    assert options.timeout == 12.5
    assert options.mode == "chembl"


def test_main__limit_zero_short_circuits(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    called = False

    def fake_run_cli_command(*_: object, **__: object) -> int:
        nonlocal called
        called = True
        return 0

    @contextmanager
    def fake_setup_cli_logging(*_: object, **__: object):
        yield SimpleNamespace(log_cfg=None)

    input_csv = tmp_path / "cellline.csv"
    input_csv.write_text("cell_chembl_id\nCHEMBL1\n", encoding="utf-8")

    monkeypatch.setattr(get_cellline_data, "setup_cli_logging", fake_setup_cli_logging)
    monkeypatch.setattr(get_cellline_data, "run_cli_command", fake_run_cli_command)
    monkeypatch.setattr(
        get_cellline_data,
        "cli",
        SimpleNamespace(prepare_io_paths=lambda *args, **kwargs: None),
    )

    exit_code = get_cellline_data.main(["--input", str(input_csv), "--limit", "0"])

    assert exit_code == 0
    assert called is False


def test_main__rejects_negative_offset(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    @contextmanager
    def fake_setup_cli_logging(*_: object, **__: object):
        yield SimpleNamespace(log_cfg=None)

    input_csv = tmp_path / "cellline.csv"
    input_csv.write_text("cell_chembl_id\nCHEMBL1\n", encoding="utf-8")

    monkeypatch.setattr(get_cellline_data, "setup_cli_logging", fake_setup_cli_logging)
    monkeypatch.setattr(get_cellline_data, "run_cli_command", lambda *args, **kwargs: 0)
    monkeypatch.setattr(
        get_cellline_data,
        "cli",
        SimpleNamespace(prepare_io_paths=lambda *args, **kwargs: None),
    )

    with pytest.raises(SystemExit):
        get_cellline_data.main(["--input", str(input_csv), "--offset", "-1"])
