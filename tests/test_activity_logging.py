from __future__ import annotations

from datetime import timezone
from pathlib import Path

import pytest

from scripts import get_activity_data


@pytest.mark.e2e
def test_get_activity_data__creates_log_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    input_path = workspace / "input.csv"
    input_path.write_text("activity_id\nCHEMBL1\n", encoding="utf-8")
    output_path = workspace / "output.csv"

    config_path = Path(__file__).resolve().parents[1] / "config" / "config.yaml"

    monkeypatch.chdir(workspace)
    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(workspace))
    monkeypatch.setattr(
        "library.cli.logging._current_date_str", lambda: "20240102_0000"
    )

    original_datetime = get_activity_data.datetime

    class _FixedDateTime(original_datetime):
        @classmethod
        def now(cls, tz=None):  # type: ignore[override]
            tzinfo = tz or timezone.utc
            return original_datetime(2024, 1, 2, 0, 0, tzinfo=tzinfo)

    monkeypatch.setattr(get_activity_data, "datetime", _FixedDateTime)
    monkeypatch.setattr("library.cli.logging.datetime", _FixedDateTime)

    def _stub_run_cli_command(
        *,
        args: object,
        parser: object,
        log_cfg: object,
        mapping: object,
        run: object,
        logger: object,
        base_parser: object | None = None,
    ) -> int:
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

    log_dir = workspace / "logs"
    assert log_dir.is_dir()

    log_path = log_dir / "get_activity_data_20240102_0000.log"
    assert log_path.is_file()

    content = log_path.read_text(encoding="utf-8")
    assert "Starting pipeline" in content
    assert "Export complete" in content
