"""Unit tests for :mod:`library.utils.data_correlation`."""

from __future__ import annotations

import pandas as pd

from library.table_quality import TableQualityProfiler
from library.utils.data_correlation import build_correlation_matrix


def test_build_correlation_matrix__matches_profiler_output(tmp_path):
    frame = pd.DataFrame(
        {
            "col_a": [1.0, 2.0, 3.0, 4.0],
            "col_b": [2.0, 4.0, 6.0, 8.0],
            "text": ["a", "b", "c", "d"],
        }
    )

    correlation = build_correlation_matrix(frame, table_name="example_table")

    profiler = TableQualityProfiler()
    profiler.consume(frame)
    _, expected_corr = profiler.build(
        table_name="example_table", destination_dir=tmp_path
    )

    pd.testing.assert_frame_equal(correlation, expected_corr)


def test_build_correlation_matrix__returns_empty_when_no_numeric():
    frame = pd.DataFrame({"letters": ["a", "b", "c"]})

    correlation = build_correlation_matrix(frame, table_name="letters")

    assert correlation.empty
    assert list(correlation.columns) == []
