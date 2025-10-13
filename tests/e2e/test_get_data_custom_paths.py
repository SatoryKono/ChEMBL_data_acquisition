from pathlib import Path

import pytest


@pytest.mark.e2e
def test_get_data_main__honours_custom_output_dir(tmp_path, monkeypatch):
    """Orchestrator should respect custom base/output directories."""

    from scripts import get_data as cli

    workspace = tmp_path / "workspace"
    output_relative = Path("artifacts/output")
    output_dir = workspace / output_relative
    output_dir.mkdir(parents=True)

    csv_one = output_dir / "output.test_20240101.csv"
    csv_two = output_dir / "output.test_20240102.csv"
    for path in (csv_one, csv_two):
        path.write_text("id\n1\n", encoding="utf-8")

    tmp_file = output_dir / "stale.tmp"
    tmp_file.write_text("junk", encoding="utf-8")

    raw_dir = output_dir / "raw"
    raw_dir.mkdir()
    (raw_dir / "temporary.csv").write_text("value\n1\n", encoding="utf-8")

    logs_dir = tmp_path / "logs"
    monkeypatch.setattr(cli, "LOGS_DIR", logs_dir)
    monkeypatch.setattr(cli, "STAGES", tuple())

    monkeypatch.chdir(tmp_path)

    exit_code = cli.main(
        [
            "--base-path",
            "workspace",
            "--output-dir",
            str(output_relative),
        ]
    )

    assert exit_code == 0
    assert csv_one.exists()
    assert csv_two.exists()
    assert not tmp_file.exists()
    assert not raw_dir.exists()

    log_files = sorted(logs_dir.glob("get_data_*.log"))
    assert len(log_files) == 1
    log_content = log_files[0].read_text(encoding="utf-8")

    resolved_output = output_dir.resolve()
    assert "Найдено 2 CSV-файлов" in log_content
    assert str(resolved_output) in log_content
    assert "[CLEANUP] Завершено: удалено" in log_content
