"""Unit tests for :mod:`library.utils.qc_report`."""

from __future__ import annotations

import pandas as pd
import pytest

import library.utils.qc_report as qc_report
from library.qa.table_quality import TableQualityProfiler as QaTableQualityProfiler
from library.table_quality import TableQualityProfiler as LegacyTableQualityProfiler
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

    assert summary.shape[1] == len(EXPECTED_COLUMNS)
    assert set(summary.columns) == set(EXPECTED_COLUMNS)

    columns = list(summary["column"])
    assert {"numeric", "text"}.issubset(columns)
    if "boolean_like" in columns:
        assert columns == ["boolean_like", "numeric", "text"]
    else:
        assert columns == ["numeric", "text"]


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


@pytest.mark.parametrize(
    "profiler_cls",
    [
        pytest.param(QaTableQualityProfiler, id="qa"),
        pytest.param(LegacyTableQualityProfiler, id="legacy"),
    ],
)
def test_build_qc_summary__accepts_prefilled_profiler(profiler_cls):
    frame = pd.DataFrame({"value": [1, 2, 3]})
    profiler = profiler_cls()
    profiler.consume(frame)

    direct = build_qc_summary(frame, table_name="demo")
    reuse = build_qc_summary(None, table_name="demo", profiler=profiler)

    pd.testing.assert_frame_equal(direct, reuse)


def test_build_qc_summary__accepts_profiler_when_optional_imports_missing(monkeypatch):
    frame = pd.DataFrame({"value": [1, 2, 3]})
    profiler = QaTableQualityProfiler()
    profiler.consume(frame)

    direct = build_qc_summary(frame, table_name="demo")

    monkeypatch.setattr(
        qc_report,
        "_TABLE_PROFILER_TYPES",
        (LegacyTableQualityProfiler,),
        raising=False,
    )

    reuse = qc_report.build_qc_summary(None, table_name="demo", profiler=profiler)

    pd.testing.assert_frame_equal(direct, reuse)


def test_build_qc_summary__reprofiles_when_profiler_invalid():
    frame = pd.DataFrame({"value": [1, 2, 3]})

    expected = build_qc_summary(frame, table_name="invalid")
    fallback = build_qc_summary(frame, table_name="invalid", profiler=object())

    pd.testing.assert_frame_equal(expected, fallback)
