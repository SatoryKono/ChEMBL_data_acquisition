"""Smoke tests for the ``get_activity_data`` CLI entry point."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from scripts import get_activity_data


def test_main_writes_expected_artifacts(tmp_path: Path, monkeypatch) -> None:
    dataset = pd.DataFrame(
        [
            {
                "activity_id": 1,
                "assay_chembl_id": "CHEMBL1",
                "target_chembl_id": "CHEMBLT",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 1.0,
                "standard_units": "NM",
            }
        ]
    )
    correlation = pd.DataFrame()
    quality = pd.DataFrame({"column": ["activity_id"]})

    def _run_activity_pipeline(**_: object):
        return dataset, correlation, quality

    monkeypatch.setattr(get_activity_data, "run_activity_pipeline", _run_activity_pipeline)

    config_path = tmp_path / "config.yaml"
    config_path.write_text("sources: {}\n", encoding="utf-8")
    output_dir = tmp_path / "out"

    exit_code = get_activity_data.main(
        [
            "--limit",
            "1",
            "--date-tag",
            "20250101",
            "--output-dir",
            str(output_dir),
            "--config",
            str(config_path),
        ]
    )

    assert exit_code == 0

    dataset_path = output_dir / "output.activity_20250101.csv"
    correlation_path = output_dir / "output.activity_20250101_data_correlation_report_table.csv"
    quality_path = output_dir / "output.activity_20250101_quality_report_table.csv"
    meta_path = dataset_path.with_name(dataset_path.name + ".meta.yaml")

    assert dataset_path.exists()
    assert correlation_path.exists()
    assert quality_path.exists()
    assert meta_path.exists()

    written = pd.read_csv(dataset_path)
    assert list(written.columns) == list(dataset.columns)
