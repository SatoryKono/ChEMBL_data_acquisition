"""Unit tests for :mod:`library.utils.qc_report`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.utils.qc_report import build_qc_summary


EXPECTED_COLUMNS = [
    "column",
    "non_null",
    "non_empty",
    "empty_pct",
    "unique_cnt",
    "unique_pct_of_non_empty",
    "pattern_cov_doi",
    "pattern_cov_issn",
    "pattern_cov_isbn",
    "pattern_cov_url",
    "pattern_cov_email",
    "bool_like_cov",
    "numeric_cov",
    "numeric_min",
    "numeric_p50",
    "numeric_p95",
    "numeric_max",
    "numeric_mean",
    "numeric_std",
    "date_cov",
    "date_min",
    "date_p50",
    "date_max",
    "text_len_min",
    "text_len_p50",
    "text_len_p95",
    "text_len_max",
    "guessed_roles",
    "top_values",
]


def test_build_qc_summary__produces_expected_structure():
    frame = pd.DataFrame(
        {
            "numeric": [1, 2, None, 2],
            "text": ["alpha", "beta", "", None],
        }
    )

    summary = build_qc_summary(frame, table_name="example_table")

    assert summary.shape == (3, len(EXPECTED_COLUMNS))
    assert list(summary.columns) == EXPECTED_COLUMNS
    assert list(summary["column"]) == ["boolean_like", "numeric", "text"]


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


def test_build_qc_summary__accepts_prefilled_profiler():
    frame = pd.DataFrame({"value": [1, 2, 3]})
    profiler = TableQualityProfiler()
    profiler.consume(frame)

    direct = build_qc_summary(frame, table_name="demo")
    reuse = build_qc_summary(None, table_name="demo", profiler=profiler)

    pd.testing.assert_frame_equal(direct, reuse)
