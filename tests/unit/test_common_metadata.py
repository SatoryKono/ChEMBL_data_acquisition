"""Tests for :mod:`library.common.metadata`."""

from __future__ import annotations

from pathlib import Path

import yaml

from config.paths import DICTIONARY_DIR

from library.common.metadata import Stats, file_sha256, write_meta_yaml


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
        config_subset={},
        inputs={},
        stats=stats,
        schema="TestSchema",
        dictionary_resources=("dictionary_root", "target_types"),
    )

    assert meta_path.exists()

    meta = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert isinstance(meta, dict)
    dictionaries = meta.get("dictionaries")
    assert isinstance(dictionaries, dict)

    manifest = yaml.safe_load((DICTIONARY_DIR / "manifest.yaml").read_text(encoding="utf-8"))
    resources = manifest.get("resources", {}) if isinstance(manifest, dict) else {}

    root_entry = resources.get("dictionary_root", {}) if isinstance(resources, dict) else {}
    target_entry = resources.get("target_types", {}) if isinstance(resources, dict) else {}

    assert dictionaries.get("dictionary_root") == {
        "version": root_entry.get("version"),
        "sha256": root_entry.get("sha256"),
    }
    assert dictionaries.get("target_types") == {
        "version": target_entry.get("version"),
        "sha256": target_entry.get("sha256"),
    }
