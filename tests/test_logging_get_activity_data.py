from __future__ import annotations

from pathlib import Path

import pytest

from scripts import get_activity_data


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
        log_path: Path | None = None,
        log_cfg: object,
        **kwargs: object,
    ) -> int:
        get_activity_data.configure_logger(log_cfg)
        if log_path is not None:
            get_activity_data.logger.info("log_destination", path=str(log_path))
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

    lines = content.splitlines()
    assert any("log_destination" in line for line in lines)
    assert any(
        "Starting get_activity_data run" in line and str(input_path) in line
        for line in lines
    )
    assert any(
        "Exported activities" in line and "rows=2" in line for line in lines
    )
    assert any(
        "Completed get_activity_data run" in line and str(output_path) in line
        for line in lines
    )
