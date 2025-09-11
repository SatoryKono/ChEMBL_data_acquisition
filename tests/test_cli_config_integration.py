from __future__ import annotations

from pathlib import Path

from scripts.table_quality_main import main


def test_table_quality_cli_with_config(tmp_path: Path) -> None:
    csv_path = tmp_path / "data.csv"
    csv_path.write_text("a\n1\n")
    config = Path("tests/fixtures/config.min.yaml")
    rc = main(
        [
            "--config",
            str(config),
            "--input",
            str(csv_path),
            "--table-name",
            "demo",
            "--output",
            str(tmp_path),
        ]
    )
    assert rc == 0
