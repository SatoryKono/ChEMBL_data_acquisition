"""Unit tests for :mod:`library.utils.data_correlation`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.qa.table_quality import TableQualityProfiler as QaTableQualityProfiler
from library.table_quality import TableQualityProfiler as LegacyTableQualityProfiler
from library.utils.data_correlation import build_correlation_matrix


@pytest.fixture()
def multi_numeric_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "col_a": [1.0, 2.0, 3.0, 4.0],
            "col_b": [2.0, 4.0, 6.0, 8.0],
            "col_c": [4.0, 3.0, 2.0, 1.0],
            "text": ["a", "b", "c", "d"],
        }
    )


def test_build_correlation_matrix__returns_expected_labels_and_shape(multi_numeric_frame):
    correlation = build_correlation_matrix(
        multi_numeric_frame, table_name="example_table"
    )

    assert correlation.shape == (3, 3)
    assert correlation.columns.tolist() == ["col_a", "col_b", "col_c"]
    assert correlation.index.tolist() == ["col_a", "col_b", "col_c"]
    assert correlation.loc["col_a", "col_c"] == pytest.approx(-1.0)
    assert correlation.loc["col_a", "col_b"] == pytest.approx(1.0)


def test_build_correlation_matrix__returns_empty_for_empty_frame():
    frame = pd.DataFrame()

    correlation = build_correlation_matrix(frame, table_name="empty")

    assert correlation.empty
    assert correlation.columns.tolist() == []


@pytest.mark.parametrize(
    "profiler_cls",
    [
        pytest.param(QaTableQualityProfiler, id="qa"),
        pytest.param(LegacyTableQualityProfiler, id="legacy"),
    ],
)
def test_build_correlation_matrix__accepts_prefilled_profiler(profiler_cls):
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [0.5, 1.0, 1.5]})
    profiler = profiler_cls()
    profiler.consume(frame)

    direct = build_correlation_matrix(frame, table_name="demo")
    reuse = build_correlation_matrix(None, table_name="demo", profiler=profiler)

    pd.testing.assert_frame_equal(direct, reuse)


def test_build_correlation_matrix__accepts_duck_typed_profiler():
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [0.5, 1.0, 1.5]})
    profiler = LegacyTableQualityProfiler()
    profiler.consume(frame)

    direct = build_correlation_matrix(frame, table_name="demo")

    class ProfilerProxy:
        def __init__(self, delegate):
            self._delegate = delegate
            self._columns = delegate._columns
            self._accumulators = delegate._accumulators

        def __getattr__(self, name):
            return getattr(self._delegate, name)

    proxy = ProfilerProxy(profiler)

    reuse = build_correlation_matrix(None, table_name="demo", profiler=proxy)

    pd.testing.assert_frame_equal(direct, reuse)


def test_build_correlation_matrix__reprofiles_when_profiler_invalid(multi_numeric_frame):
    expected = build_correlation_matrix(multi_numeric_frame, table_name="fallback")

    fallback = build_correlation_matrix(
        multi_numeric_frame,
        table_name="fallback",
        profiler=object(),
    )

    pd.testing.assert_frame_equal(expected, fallback)
def test_build_correlation_matrix__falls_back_when_profiler_invalid(
    multi_numeric_frame,
):
    class DummyProfiler:
        pass

    with pytest.warns(RuntimeWarning, match="Ignoring incompatible profiler"):
        correlation = build_correlation_matrix(
            multi_numeric_frame,
            table_name="fallback",
            profiler=DummyProfiler(),
        )

    assert correlation.shape == (3, 3)
