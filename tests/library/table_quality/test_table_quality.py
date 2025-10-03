from __future__ import annotations

import warnings
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from library.config import Config
from library.qa.table_quality import TableQualityProfiler, analyze_table_quality
from library.utils.cli_tools.table_quality_main import run


def test_analyze_table_quality(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "num": [1, 2, 3, 4],
            "str": ["10.1234/abc", "", None, "test"],
            "flag": pd.Series([True, False, True, False], dtype="boolean"),
        }
    )
    quality, corr = analyze_table_quality(
        df,
        table_name="sample",
        destination_dir=tmp_path,
    )

    assert set(quality["column"]) == {"num", "str", "flag"}
    assert (tmp_path / "sample_quality_report_table.csv").exists()
    assert (tmp_path / "sample_data_correlation_report_table.csv").exists()
    assert corr.shape == (2, 2)


def test_profiler_build_sanitizes_table_name(tmp_path: Path) -> None:
    profiler = TableQualityProfiler()
    profiler.consume(pd.DataFrame({"value": [1, 2]}))

    nested_name = "nested/dir/output.csv"
    profiler.build(table_name=nested_name, destination_dir=tmp_path)

    expected_quality = tmp_path / "output.csv_quality_report_table.csv"
    expected_corr = tmp_path / "output.csv_data_correlation_report_table.csv"

    assert expected_quality.exists()
    assert expected_corr.exists()


def test_analyze_table_quality_supports_relative_destination(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    df = pd.DataFrame({"value": [1, 2, 3]})
    monkeypatch.chdir(tmp_path)
    destination = Path("relative_output")

    quality, _ = analyze_table_quality(
        df,
        table_name="relative",
        destination_dir=destination,
    )

    expected_dir = tmp_path / destination
    assert expected_dir.exists()
    assert (expected_dir / "relative_quality_report_table.csv").exists()
    assert not quality.empty


def test_analyze_table_quality_suppresses_warnings(tmp_path: Path) -> None:
    csv_path = tmp_path / "mixed.csv"
    csv_path.write_text("num,str\n1,2020-01-01\nx,1a; DPCPX\n")

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        analyze_table_quality(
            csv_path,
            table_name="mixed",
            destination_dir=tmp_path,
        )

    assert (tmp_path / "mixed_quality_report_table.csv").exists()


def test_analyze_table_quality_handles_sequences(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "seq": [[1], [], np.array([1]), np.array([])],
        }
    )
    quality, _ = analyze_table_quality(
        df,
        table_name="seq",
        destination_dir=tmp_path,
    )

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


def test_table_quality_run_handles_mixed_types(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "mixed": [1, "2", 3.5, None],
            "identifier": ["0001", "0002", "ABC123", None],
        }
    )
    csv_path = tmp_path / "mixed.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=tmp_path,
        table_name="mixed",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)
    assert exit_code == 0

    report_path = tmp_path / "mixed_quality_report_table.csv"
    assert report_path.exists()
    report = pd.read_csv(report_path)

    expected_quality, _ = analyze_table_quality(
        df,
        table_name="expected",
        destination_dir=tmp_path,
    )

    expected_mixed = expected_quality.loc[expected_quality["column"] == "mixed"].iloc[0]
    actual_mixed = report.loc[report["column"] == "mixed"].iloc[0]

    assert actual_mixed["numeric_cov"] == pytest.approx(expected_mixed["numeric_cov"])
    assert actual_mixed["numeric_mean"] == pytest.approx(expected_mixed["numeric_mean"])


def test_table_quality_handles_duplicate_columns(tmp_path: Path) -> None:
    profiler = TableQualityProfiler()
    df = pd.DataFrame([[1, 10], [None, 20]], columns=["dup", "dup"])

    profiler.consume(df)
    quality, _ = profiler.build(table_name="duplicates", destination_dir=tmp_path)

    dup_row = quality.loc[quality["column"] == "dup"].iloc[0]
    assert int(dup_row["non_null"]) == 1
    assert pytest.approx(float(dup_row["empty_pct"])) == 0.5


def test_table_quality_run_rejects_output_with_suffix(tmp_path: Path) -> None:
    df = pd.DataFrame({"value": [1, 2, 3]})
    csv_path = tmp_path / "input.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()

    output_path = tmp_path / "report.csv"

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=output_path,
        table_name="invalid",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)

    assert exit_code == 1
    assert not output_path.exists()


def test_table_quality_run_requires_existing_directory_when_disallowed(
    tmp_path: Path,
) -> None:
    df = pd.DataFrame({"value": [1, 2, 3]})
    csv_path = tmp_path / "input.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()
    cfg.io.exist_ok = False

    output_dir = tmp_path / "reports"

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=output_dir,
        table_name="deny",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)

    assert exit_code == 1
    assert not output_dir.exists()


def test_table_quality_run_creates_directory_when_allowed(tmp_path: Path) -> None:
    df = pd.DataFrame({"value": [1, 2, 3]})
    csv_path = tmp_path / "input.csv"
    df.to_csv(csv_path, index=False)

    cfg = Config()
    cfg.io.exist_ok = True

    output_dir = tmp_path / "reports"

    args = Namespace(
        input_csv=csv_path,
        sep=",",
        encoding="utf-8-sig",
        output_csv=output_dir,
        table_name="created",
        doc_quality_enable=None,
        sample_rows=None,
        include_columns=None,
        exclude_columns=None,
    )

    exit_code = run(cfg, args)

    assert exit_code == 0
    assert output_dir.is_dir()
    report = output_dir / "created_quality_report_table.csv"
    correlation = output_dir / "created_data_correlation_report_table.csv"
    assert report.exists()
    assert correlation.exists()
