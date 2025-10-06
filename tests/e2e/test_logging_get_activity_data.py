from __future__ import annotations

from pathlib import Path

import importlib
import importlib.metadata as metadata

import pytest


def _load_cli_module(entry_point: str):
    entries = metadata.entry_points(group="console_scripts", name=entry_point)
    if not entries:
        msg = f"console script '{entry_point}' is not registered"
        raise LookupError(msg)
    module_path, _, attribute = entries[0].value.partition(":")
    if module_path.startswith("library.cli.commands."):
        from library.cli.commands import _resolve_module

        module = _resolve_module(module_path.rsplit(".", 1)[-1])
    else:
        module = importlib.import_module(module_path)
    if attribute and not hasattr(module, attribute):  # pragma: no cover - guard rail
        msg = f"entry point '{entry_point}' refers to missing attribute '{attribute}'"
        raise AttributeError(msg)
    return module


get_activity_data = _load_cli_module("get-activity-data")


@pytest.mark.e2e
def test_logging_get_activity_data__writes_expected_messages(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_dir = tmp_path / "cli"
    base_dir.mkdir()

    input_path = base_dir / "input.csv"
    input_path.write_text("activity_id\nCHEMBL1\n", encoding="utf-8")
    input_path = input_path.resolve()

    output_path = (base_dir / "output.csv").resolve()

    config_path = base_dir / "config.yaml"
    config_path.write_text(
        "io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n",
        encoding="utf-8",
    )

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_dir))
    monkeypatch.setattr(
        "library.cli.logging._current_date_str", lambda: "20240102"
    )

    class _FixedUUID:
        hex = "run-id-0001"

    monkeypatch.setattr(
        "library.cli.parser.uuid.uuid4", lambda: _FixedUUID()
    )

    def _stub_run_cli_command(
        *,
        args: object,
        parser: object,
        log_cfg: object,
        **kwargs: object,
    ) -> int:
        get_activity_data.configure_logger(log_cfg)
        get_activity_data.logger.info(
            "Starting get_activity_data run",
            input=str(getattr(args, "input_csv", "")),
            output=str(getattr(args, "final_out", "")),
        )
        get_activity_data.logger.info(
            "Exported activities",
            output=str(getattr(args, "final_out", "")),
            rows=2,
        )
        get_activity_data.logger.info(
            "Completed get_activity_data run",
            output=str(getattr(args, "final_out", "")),
        )
        return 0

    monkeypatch.setattr(get_activity_data, "run_cli_command", _stub_run_cli_command)
    monkeypatch.setattr("library.cli_utils.run_cli_command", _stub_run_cli_command)

    exit_code = get_activity_data.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--final-out",
            str(output_path),
        ]
    )
    assert exit_code == 0

    log_dir = base_dir / "logs"
    assert log_dir.is_dir()

    log_files = sorted(log_dir.glob("*.log"))
    assert len(log_files) == 1

    log_path = log_files[0]
    assert log_path.name == "get_activity_data_20240102.log"

    content = log_path.read_text(encoding="utf-8")

    assert "Starting get_activity_data run" in content
    assert "Exported activities" in content
    assert "Completed get_activity_data run" in content

    expected_lines = [
        (
            "[2020-01-01 00:00:00,000] [INFO] [chembl] "
            f"Starting get_activity_data run input='{input_path}' output='{output_path}' "
            "rps=None run_id='run-id-0001' status=None"
        ),
        (
            "[2020-01-01 00:00:00,000] [INFO] [chembl] "
            f"Exported activities output='{output_path}' rows=2 "
            "rps=None run_id='run-id-0001' status=None"
        ),
        (
            "[2020-01-01 00:00:00,000] [INFO] [chembl] "
            f"Completed get_activity_data run output='{output_path}' "
            "rps=None run_id='run-id-0001' status=None"
        ),
    ]

    assert content.splitlines() == expected_lines
