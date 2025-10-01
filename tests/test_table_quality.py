from __future__ import annotations

import os
import warnings
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd

from library.config import Config
from library.table_quality import analyze_table_quality
from library.utils.cli_tools.table_quality_main import run


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


def test_analyze_table_quality_handles_sequences(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "seq": [[1], [], np.array([1]), np.array([])],
        }
    )
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        quality, _ = analyze_table_quality(df, table_name="seq")
    finally:
        os.chdir(cwd)

    non_empty = int(quality.loc[quality["column"] == "seq", "non_empty"].iloc[0])
    assert non_empty == 2


def test_table_quality_run_skips_when_disabled(tmp_path: Path) -> None:
    df = pd.DataFrame({"keep": [1, 2], "drop": [3, 4]})
    csv_path = tmp_path / "input.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()
    cfg.system.doc_quality.enable = False

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=tmp_path,
        table_name="disabled",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)

    assert exit_code == 0
    assert not any(tmp_path.glob("disabled_*"))


def test_table_quality_run_respects_sampling_and_filters(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "keep": [1, 2, 3, 4],
            "drop": [10, 20, 30, 40],
            "text": ["a", "b", "c", "d"],
        }
    )
    csv_path = tmp_path / "filtered.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()
    cfg.system.doc_quality.sample_rows = 2
    cfg.system.doc_quality.include_columns = ("keep", "drop")
    cfg.system.doc_quality.exclude_columns = ("drop",)

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=tmp_path,
        table_name="filtered",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)
    assert exit_code == 0

    report_path = tmp_path / "filtered_quality_report_table.csv"
    assert report_path.exists()

    report = pd.read_csv(report_path)
    assert set(report["column"]) == {"keep"}
    non_null = int(report.loc[report["column"] == "keep", "non_null"].iloc[0])
    assert non_null == 2
