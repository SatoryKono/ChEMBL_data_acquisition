from __future__ import annotations

import yaml
from pathlib import Path

from library.metadata import write_meta_yaml, Stats


def test_write_meta_yaml_creates_file(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.csv"
    csv_path.write_text("a,b\n1,2\n", encoding="utf-8")

    stats: Stats = {
        "rows_total": 2,
        "rows_kept": 2,
        "rows_dropped": 0,
        "output_sha256": "deadbeef",
    }

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="unit-test",
        config_subset={"api_key": "secret"},
        inputs={"source": "dummy"},
        stats=stats,
        schema="TestSchema",
    )

    assert meta_path.exists()
    with meta_path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    assert data["command"] == "unit-test"
    assert data["config"]["api_key"] == "***"
    assert data["inputs"]["source"] == "dummy"
    assert data["stats"] == stats
    assert data["schema"] == "TestSchema"
    assert "generated_at" in data
    assert "git_sha" in data
