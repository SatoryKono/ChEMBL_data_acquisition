"""Tests for :mod:`library.io.metadata`."""

from __future__ import annotations

from pathlib import Path

import yaml

from library.io.metadata import write_meta_yaml
from library.pipelines.common.metadata import get_pipeline_version


def test_write_meta_yaml__includes_pipeline_version(tmp_path: Path) -> None:
    csv_path = tmp_path / "output.csv"
    csv_path.write_text("value\n1\n", encoding="utf-8")

    meta_path = write_meta_yaml(csv_path, columns=["value"], dtypes={"value": "Int64"})

    assert meta_path.exists()

    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert isinstance(meta, dict)
    assert meta.get("pipeline_version") == get_pipeline_version()
