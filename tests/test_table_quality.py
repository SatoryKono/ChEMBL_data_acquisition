from __future__ import annotations

import os
import warnings
from pathlib import Path

import pandas as pd

from library.table_quality import analyze_table_quality


def test_analyze_table_quality(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "num": [1, 2, 3, 4],
            "str": ["10.1234/abc", "", None, "test"],
            "flag": pd.Series([True, False, True, False], dtype="boolean"),
        }
    )
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        quality, corr = analyze_table_quality(df, table_name="sample")
    finally:
        os.chdir(cwd)

    assert set(quality["column"]) == {"num", "str", "flag"}
    assert (tmp_path / "sample_quality_report_table.csv").exists()
    assert (tmp_path / "sample_data_correlation_report_table.csv").exists()
    assert corr.shape == (2, 2)


def test_analyze_table_quality_suppresses_warnings(tmp_path: Path) -> None:
    csv_path = tmp_path / "mixed.csv"
    csv_path.write_text("num,str\n1,2020-01-01\nx,1a; DPCPX\n")

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            analyze_table_quality(csv_path, table_name="mixed")
    finally:
        os.chdir(cwd)

    assert (tmp_path / "mixed_quality_report_table.csv").exists()
