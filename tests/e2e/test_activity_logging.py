from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from scripts import get_activity_data


@pytest.mark.e2e
def test_activity_logging__relative_env_base_anchored_to_repo_root(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    repo_root = Path(__file__).resolve().parents[2]
    scripts_dir = repo_root / "scripts"
    monkeypatch.chdir(scripts_dir)

    relative_base = "activity_logs_tmp"
    target_base = repo_root / relative_base
    log_dir = target_base / "logs"
    if target_base.exists():
        shutil.rmtree(target_base)

    input_path = tmp_path / "input.csv"
    input_path.write_text("activity_id\nCHEMBL1\n", encoding="utf-8")

    output_path = tmp_path / "output.csv"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n", encoding="utf-8"
    )

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", relative_base)
    monkeypatch.setattr("library.cli.logging._current_date_str", lambda: "20240102")

    def _stub_run_cli_command(
        *,
        args: object,
        parser: object,
        log_cfg: object,
        **_: object,
    ) -> int:
        get_activity_data.configure_logger(log_cfg)
        get_activity_data.logger.info(
            "activity_pipeline_start",
            input=str(getattr(args, "input_csv", "")),
        )
        get_activity_data.logger.info(
            "activity_pipeline_done",
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

    try:
        assert exit_code == 0
        log_path = log_dir / "get_activity_data_20240102.log"
        assert log_path.exists()
        content = log_path.read_text(encoding="utf-8")
        assert "activity_pipeline_start" in content
        assert "activity_pipeline_done" in content
    finally:
        if target_base.exists():
            shutil.rmtree(target_base)
