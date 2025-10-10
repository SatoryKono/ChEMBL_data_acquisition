"""Tests for :mod:`library.common.metadata`."""

from __future__ import annotations

from pathlib import Path

import yaml

from library.common.metadata import Stats, file_sha256, write_meta_yaml
from library.pipelines.common.metadata import get_pipeline_version
from library.resources.dictionaries import get_resource


def test_write_meta_yaml__records_dictionary_metadata(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.csv"
    csv_path.write_text("value\n1\n", encoding="utf-8")

    stats: Stats = {
        "rows_total": 1,
        "rows_kept": 1,
        "rows_dropped": 0,
        "output_sha256": file_sha256(csv_path),
    }

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command="pytest",
        config={},
        inputs={},
        stats=stats,
        schema="TestSchema",
        dictionary_resources=("dictionary_root", "target_types"),
    )

    assert meta_path.exists()

    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert isinstance(meta, dict)
    assert meta.get("pipeline_version") == get_pipeline_version()
    dictionaries = meta.get("dictionaries")
    assert isinstance(dictionaries, dict)

    root_resource = get_resource("dictionary_root")
    target_resource = get_resource("target_types")

    assert dictionaries.get("dictionary_root") == {
        "version": root_resource.version,
        "sha256": root_resource.sha256,
    }
    assert dictionaries.get("target_types") == {
        "version": target_resource.version,
        "sha256": target_resource.sha256,
    }
