"""End-to-end checks for CLI logging helpers."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pytest

from library.cli.base import PipelineCLIBase
from scripts import (
    get_activity_data,
    get_assay_data,
    get_document_data,
    get_target_data,
    get_testitem_data,
)
from tests.helpers.logs import parse_log_file


def _program_name_from_module(module: Any) -> str:
    cli_candidate = getattr(module, "_CLI", None)
    if isinstance(cli_candidate, PipelineCLIBase):
        return cli_candidate.get_program_name()
    program_name = getattr(module, "PROGRAM_NAME", None)
    if isinstance(program_name, str) and program_name:
        return program_name
    module_file = getattr(module, "__file__", None)
    if module_file:
        return Path(module_file).with_suffix("").name
    return module.__name__.rsplit(".", 1)[-1]


_SCRIPT_CASES = (
    {
        "module": get_activity_data,
        "command": None,
        "input_flag": "--input",
        "output_flag": "--final-out",
        "output_attr": "final_out",
        "prefix": "activity_pipeline",
        "extra_args": (),
    },
    {
        "module": get_assay_data,
        "command": None,
        "input_flag": "--input",
        "output_flag": "--final-out",
        "output_attr": "final_out",
        "prefix": "assay_pipeline",
        "extra_args": (),
    },
    {
        "module": get_document_data,
        "command": "pubmed",
        "input_flag": "--input",
        "output_flag": "--final-out",
        "output_attr": "final_out",
        "prefix": "document_pipeline",
        "extra_args": (),
    },
    {
        "module": get_target_data,
        "command": "chembl",
        "input_flag": "--input",
        "output_flag": "--final-out",
        "output_attr": "final_out",
        "prefix": "target_pipeline",
        "extra_args": ("--date", "20240102"),
    },
    {
        "module": get_testitem_data,
        "command": None,
        "input_flag": "--input",
        "output_flag": "--final-out",
        "output_attr": "final_out",
        "prefix": "testitem_pipeline",
        "extra_args": (),
    },
)


def _run_logging_case(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    case: dict[str, Any],
    *,
    expected_date: str,
    date_override: str | None,
    datetime_cls: type[datetime] | None = None,
) -> tuple[Path, list[dict[str, Any]]]:
    base_path = tmp_path
    log_dir = base_path / "logs"
    log_dir.mkdir(parents=True)
    monkeypatch.setattr("library.cli.logging._DEFAULT_LOG_DIR", log_dir)

    config_path = base_path / "config.yaml"
    config_path.write_text(
        "io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n", encoding="utf-8"
    )

    input_path = base_path / "input.csv"
    input_path.write_text("identifier\nCHEMBL1\n", encoding="utf-8")

    output_path = base_path / "output.csv"

    monkeypatch.chdir(base_path)
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))
    if date_override is not None:
        monkeypatch.setattr(
            "library.cli.logging._current_date_str", lambda: date_override
        )
    elif datetime_cls is not None:
        monkeypatch.setattr("library.cli.logging.datetime", datetime_cls, raising=False)

    module = case["module"]
    prefix = case["prefix"]
    command = case["command"]
    input_flag = case["input_flag"]
    output_flag = case["output_flag"]
    output_attr = case["output_attr"]
    extra_args = list(case["extra_args"])

    if datetime_cls is not None:
        monkeypatch.setattr(module, "datetime", datetime_cls, raising=False)

    def _run_cli_command_stub(
        *,
        args: Any,
        parser: Any,
        log_cfg: Any,
        mapping: Any,
        run: Any,
        logger: Any,
        base_parser: Any | None = None,
    ) -> int:
        module.configure_logger(log_cfg)
        logger.info(
            "pipeline_start",
            run_id=log_cfg.run_id,
            command=getattr(args, "command", command or Path(module.__file__).stem),
        )
        logger.info(
            "pipeline_parameters",
            input=str(getattr(args, "input_csv", "")),
            output=str(getattr(args, output_attr, "")),
            limit=getattr(args, "limit", None),
            offset=getattr(args, "offset", None),
        )
        exit_code = run(object(), args)
        if exit_code == 0:
            logger.info("pipeline_done", run_id=log_cfg.run_id)
        else:
            logger.info("pipeline_fail", run_id=log_cfg.run_id, exit_code=exit_code)
        return int(exit_code)

    monkeypatch.setattr(module, "run_cli_command", _run_cli_command_stub, raising=False)

    def _stub_run(_cfg: Any, args: Any) -> int:
        module.logger.info(
            f"{prefix}_start",
            stage="ingest",
            input=str(getattr(args, "input_csv", "")),
        )
        module.logger.info(
            f"{prefix}_records",
            processed=2,
            discarded=1,
            stage="ingest",
        )
        module.logger.info(
            f"{prefix}_done",
            output=str(getattr(args, output_attr, "")),
        )
        return 0

    monkeypatch.setattr(module, "run", _stub_run, raising=False)

    argv: list[str] = []
    if command:
        argv.append(command)
    argv.extend(["--config", str(config_path)])
    if input_flag:
        argv.extend([input_flag, str(input_path)])
    if output_flag:
        target_output = (
            output_path if output_flag != "--final-out" else base_path / "targets.csv"
        )
        argv.extend([output_flag, str(target_output)])
    argv.extend(extra_args)

    exit_code = module.main(argv)
    assert exit_code == 0

    program_name = _program_name_from_module(module)
    expected_prefix = program_name
    log_files = sorted(
        path for path in log_dir.glob("*.log") if path.name.startswith(expected_prefix)
    )
    if not log_files:
        fallback_dir = base_path / "logs"
        log_files = sorted(
            path
            for path in fallback_dir.glob("*.log")
            if path.name.startswith(expected_prefix)
        )
    assert len(log_files) == 1
    log_path = log_files[0]
    expected_name = f"{program_name}_{expected_date}.log"
    assert log_path.name == expected_name

    events = parse_log_file(log_path)
    return log_path, events


@pytest.mark.e2e
@pytest.mark.parametrize("case", _SCRIPT_CASES, ids=lambda c: _program_name_from_module(c["module"]))
def test_cli_logging__creates_log_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, case: dict[str, Any]
) -> None:
    class _FixedDateTime(datetime):
        @classmethod
        def now(cls, tz=None):  # type: ignore[override]
            tzinfo = tz or UTC
            return datetime(2024, 1, 2, 0, 0, tzinfo=tzinfo)

    datetime_cls = _FixedDateTime if case["module"] is get_activity_data else None

    _, events = _run_logging_case(
        tmp_path,
        monkeypatch,
        case,
        expected_date="20240102",
        date_override="20240102",
        datetime_cls=datetime_cls,
    )

    prefix = case["prefix"]
    event_names = {record.get("event") for record in events}
    expected_events = {
        "pipeline_start",
        "pipeline_parameters",
        f"{prefix}_start",
        f"{prefix}_records",
        f"{prefix}_done",
        "pipeline_done",
    }
    assert expected_events.issubset(event_names)

    record_entry = next(
        record for record in events if record.get("event") == f"{prefix}_records"
    )
    data = record_entry.get("data", {})
    assert data.get("processed") == 2
    assert data.get("discarded") == 1


@pytest.mark.e2e
def test_cli_logging__uses_datetime_hook_for_default_date(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = next(item for item in _SCRIPT_CASES if item["module"] is get_activity_data)

    class _FixedDateTime(datetime):
        @classmethod
        def now(cls, tz=None):  # type: ignore[override]
            tzinfo = tz or UTC
            return datetime(2024, 1, 2, 0, 0, tzinfo=tzinfo)

    expected_date = _FixedDateTime.now(UTC).strftime("%Y%m%d")

    log_path, events = _run_logging_case(
        tmp_path,
        monkeypatch,
        case,
        expected_date=expected_date,
        date_override=None,
        datetime_cls=_FixedDateTime,
    )

    assert log_path.name.endswith(f"_{expected_date}.log")

    event_names = {record.get("event") for record in events}
    assert "pipeline_start" in event_names
    assert "pipeline_done" in event_names
