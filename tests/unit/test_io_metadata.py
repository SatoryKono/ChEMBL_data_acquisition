"""Tests for :mod:`library.io.metadata`."""

from __future__ import annotations

from pathlib import Path

import yaml

from library.common.metadata import Stats, file_sha256, write_meta_yaml as write_common_meta
from library.io.metadata import write_meta_yaml
from library.pipelines.common.metadata import get_pipeline_version


def _make_stats(csv_path: Path) -> Stats:
    return {
        "rows_total": 1,
        "rows_kept": 1,
        "rows_dropped": 0,
        "output_sha256": file_sha256(csv_path),
    }


def test_write_meta_yaml__includes_pipeline_version(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.csv"
    csv_path.write_text("value\n1\n", encoding="utf-8")

    meta_path = write_meta_yaml(csv_path, columns=["value"], dtypes={"value": "Int64"})

    assert meta_path.exists()

    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert isinstance(meta, dict)
    assert meta.get("pipeline_version") == get_pipeline_version()


def test_write_meta_yaml__aligns_with_common_writer(tmp_path: Path) -> None:
    csv_simple = tmp_path / "simple.csv"
    csv_simple.write_text("value\n1\n", encoding="utf-8")

    csv_detailed = tmp_path / "detailed.csv"
    csv_detailed.write_text("value\n1\n", encoding="utf-8")

    simple_meta_path = write_meta_yaml(
        csv_simple,
        columns=["value"],
        dtypes={"value": "Int64"},
    )
    detailed_meta_path = write_common_meta(
        csv_path=csv_detailed,
        command="pytest",
        config_subset={},
        inputs={"source": "unit"},
        stats=_make_stats(csv_detailed),
        schema="TestSchema",
    )

    simple_meta = yaml.safe_load(simple_meta_path.read_text(encoding="utf-8"))
    detailed_meta = yaml.safe_load(detailed_meta_path.read_text(encoding="utf-8"))

    assert isinstance(simple_meta, dict)
    assert isinstance(detailed_meta, dict)

    assert set(simple_meta) == set(detailed_meta)
    assert simple_meta["schema"] is None
    assert simple_meta["stats"] == {}
    assert detailed_meta["schema"] == "TestSchema"
    assert detailed_meta["stats"]["rows_total"] == 1
