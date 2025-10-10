"""Unit tests for :mod:`library.utils.qc_report`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.utils.qc_report import build_qc_summary


def test_build_qc_summary__computes_expected_metrics():
    frame = pd.DataFrame(
        {
            "numeric": [1, 2, None, 2],
            "text": ["alpha", "beta", "", None],
        }
    )

    summary = build_qc_summary(frame, table_name="example_table")

    expected_columns = [
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
    assert summary.columns.tolist() == expected_columns

    numeric_row = summary.loc[summary["column"] == "numeric"].iloc[0]
    assert numeric_row["non_null"] == 3
    assert numeric_row["non_empty"] == 3
    assert numeric_row["empty_pct"] == pytest.approx(0.25)
    assert numeric_row["unique_cnt"] == 2
    assert numeric_row["numeric_cov"] == pytest.approx(0.75)
    assert numeric_row["guessed_roles"] == "free-text"

    text_row = summary.loc[summary["column"] == "text"].iloc[0]
    assert text_row["non_null"] == 3
    assert text_row["non_empty"] == 2
    assert text_row["empty_pct"] == pytest.approx(0.5)
    assert text_row["unique_cnt"] == 3
    assert text_row["guessed_roles"] == "identifier-like"
    assert text_row["top_values"] == "alpha (1); beta (1);  (1)"


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
