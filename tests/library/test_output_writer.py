from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.io import output_writer


def test_save_standard_outputs_orders_columns_and_prunes_sidecars(tmp_path: Path) -> None:
    dataset = pd.DataFrame(
        {
            "identifier": ["row-2", "row-1"],
            "secondary": ["B", "A"],
            "nullable": ["", None],
        }
    )
    correlation = pd.DataFrame({"metric": ["pearson"], "value": ["0.9"]})
    quality = pd.DataFrame({"column": ["identifier"], "non_null": ["2"]})

    working_dir = tmp_path / "exports"
    working_dir.mkdir(parents=True, exist_ok=True)

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="sample",
        date_tag="20240101",
        column_order=("identifier", "secondary", "nullable"),
        output_dir=working_dir,
    )

    expected_dataset = working_dir / "output.sample_20240101.csv"
    expected_corr = (
        working_dir
        / "output.sample_20240101_data_correlation_report_table.csv"
    )
    expected_quality = (
        working_dir / "output.sample_20240101_quality_report_table.csv"
    )

    assert artifacts.dataset == expected_dataset
    assert artifacts.correlation_report == expected_corr
    assert artifacts.quality_report == expected_quality

    written = pd.read_csv(expected_dataset, dtype="string")
    assert list(written.columns) == ["identifier", "secondary", "nullable"]
    assert list(written["identifier"]) == ["row-1", "row-2"]
    assert pd.isna(written.loc[1, "nullable"])

    for produced in (expected_dataset, expected_corr, expected_quality):
        assert produced.exists()
        assert not produced.with_suffix(produced.suffix + ".meta.yaml").exists()


def test_save_standard_outputs_respects_cleanup_flag(tmp_path: Path) -> None:
    dataset = pd.DataFrame({"identifier": ["row-1"]})
    correlation = pd.DataFrame({"metric": ["pearson"], "value": ["1.0"]})
    quality = pd.DataFrame({"column": ["identifier"], "non_null": ["1"]})

    legacy_dir = tmp_path / "legacy"
    legacy_dir.mkdir(parents=True, exist_ok=True)
    legacy_path = legacy_dir / "legacy-output.csv"
    legacy_path.write_text("identifier\nlegacy\n", encoding="utf-8")

    artifacts = output_writer.save_standard_outputs(
        dataset,
        correlation,
        quality,
        table_name="legacy",
        date_tag="20240202",
        output_path=legacy_path,
        cleanup_source=False,
    )

    assert legacy_path.exists(), "legacy file should remain when cleanup disabled"

    expected_dataset = legacy_dir / "output.legacy_20240202.csv"
    assert artifacts.dataset == expected_dataset
    assert expected_dataset.exists()
