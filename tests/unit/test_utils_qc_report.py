"""Unit tests for :mod:`library.utils.qc_report`."""

from __future__ import annotations

import pandas as pd

from library.table_quality import TableQualityProfiler
from library.utils.qc_report import build_qc_summary


def test_build_qc_summary__matches_profiler_output(tmp_path):
    frame = pd.DataFrame(
        {
            "numeric": [1, 2, 3, 4],
            "text": ["alpha", "beta", "", None],
            "boolean_like": ["true", "false", "true", "no"],
        }
    )

    summary = build_qc_summary(frame, table_name="example_table")

    profiler = TableQualityProfiler()
    profiler.consume(frame)
    expected_summary, _ = profiler.build(
        table_name="example_table", destination_dir=tmp_path
    )

    pd.testing.assert_frame_equal(summary, expected_summary)


def test_build_qc_summary__respects_sampling_and_filters():
    frame = pd.DataFrame(
        {
            "kept": ["a", "b", "c"],
            "ignored": [1, 2, 3],
        }
    )

    summary = build_qc_summary(
        frame,
        table_name="demo",
        include_columns=["kept"],
        sample_rows=2,
    )

    assert list(summary["column"]) == ["kept"]
    assert summary.loc[0, "non_null"] == 2
    roles_value = summary.loc[0, "guessed_roles"]
    roles = roles_value.split("|") if roles_value else []
    assert roles == ["identifier-like"]
