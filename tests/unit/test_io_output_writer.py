"""Tests for :mod:`library.io.output_writer`."""

from __future__ import annotations

import os
import random
from pathlib import Path

import numpy as np
import pandas as pd

from library.config import IoCfg
from library.io.output_writer import save_standard_outputs


def _fix_seed(seed: int = 42) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def test_save_standard_outputs__writes_expected_files(tmp_path: Path) -> None:
    _fix_seed()

    cfg = IoCfg(
        output_dir=tmp_path,
        csv_encoding="utf-8",
        csv_sep="|",
        exist_ok=True,
    )

    dataset = pd.DataFrame(
        {
            "chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
            "name": pd.Series(["Aspirin", "Ibuprofen"], dtype="string"),
        }
    )
    quality_report = pd.DataFrame(
        {
            "metric": pd.Series(["rows"], dtype="string"),
            "value": pd.Series([2], dtype="Int64"),
        }
    )
    correlation_report = pd.DataFrame(
        {
            "feature": pd.Series(["chembl_id"], dtype="string"),
            "correlated_with": pd.Series(["name"], dtype="string"),
            "score": pd.Series([1.0], dtype="Float64"),
        }
    )

    table_name = "test_table"
    date_tag = "20240101"

    artifacts = save_standard_outputs(
        dataset,
        quality_report,
        correlation_report,
        table_name=table_name,
        date_tag=date_tag,
        cfg=cfg,
    )

    stem = f"output.{table_name}_{date_tag}"
    assert artifacts.dataset == tmp_path / f"{stem}.csv"
    assert artifacts.quality_report == tmp_path / f"{stem}_quality_report_table.csv"
    assert (
        artifacts.correlation_report
        == tmp_path / f"{stem}_data_correlation_report_table.csv"
    )

    for path in (artifacts.dataset, artifacts.quality_report, artifacts.correlation_report):
        assert path.exists()

    roundtrip_dataset = pd.read_csv(
        artifacts.dataset,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        dtype="string",
    ).convert_dtypes()
    roundtrip_quality = pd.read_csv(
        artifacts.quality_report,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        dtype="string",
    ).convert_dtypes()
    roundtrip_quality = roundtrip_quality.astype({"value": "Int64"})
    roundtrip_correlation = pd.read_csv(
        artifacts.correlation_report,
        sep=cfg.csv_sep,
        encoding=cfg.csv_encoding,
        dtype="string",
    ).convert_dtypes()
    roundtrip_correlation = roundtrip_correlation.astype({"score": "Float64"})

    pd.testing.assert_frame_equal(roundtrip_dataset, dataset.convert_dtypes())
    pd.testing.assert_frame_equal(roundtrip_quality, quality_report.convert_dtypes())
    pd.testing.assert_frame_equal(
        roundtrip_correlation, correlation_report.convert_dtypes()
    )

    entries = sorted(tmp_path.iterdir())
    expected_csvs = sorted(
        [artifacts.dataset, artifacts.quality_report, artifacts.correlation_report]
    )
    expected_entries = sorted(
        expected_csvs
        + [Path(f"{path}.meta.yaml") for path in expected_csvs]
    )

    assert entries == expected_entries
    assert [path for path in entries if path.suffix == ".csv"] == expected_csvs
