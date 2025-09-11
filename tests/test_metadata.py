from __future__ import annotations

from pathlib import Path

import shutil
import yaml

from library.metadata import Stats, write_meta_yaml


def test_write_meta_yaml_creates_file(tmp_path: Path) -> None:
    data_src = Path(__file__).parent / "data" / "meta_input.csv"
    csv_path = tmp_path / "output.csv"
    shutil.copy(data_src, csv_path)

    stats: Stats = {
        "rows_total": 2,
        "rows_kept": 2,
        "rows_dropped": 0,
        "output_sha256": "deadbeef",
    }

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="unit-test",
        config_subset={"api_key": "secret", "password": "topsecret"},
        inputs={"source": "dummy"},
        stats=stats,
        schema="TestSchema",
    )

    assert meta_path.exists()
    with meta_path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    assert data["command"] == "unit-test"
    assert data["config"]["api_key"] == "***"
    assert data["config"]["password"] == "***"
    assert data["inputs"]["source"] == "dummy"
    assert data["stats"] == stats
    assert data["schema"] == "TestSchema"

    required_keys = {
        "generated_at",
        "git_sha",
        "python_version",
        "platform",
        "command",
        "config",
        "inputs",
        "stats",
        "schema",
    }
    assert required_keys <= data.keys()
