"""Unit tests for :mod:`library.utils.data_correlation`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.utils.data_correlation import build_correlation_matrix


def test_build_correlation_matrix__produces_expected_structure():
    frame = pd.DataFrame(
        {
            "col_a": [1.0, 2.0, 3.0, 4.0],
            "col_b": [2.0, 4.0, 6.0, 8.0],
            "text": ["a", "b", "c", "d"],
        }
    )

    correlation = build_correlation_matrix(frame, table_name="example_table")

    assert correlation.shape == (2, 2)
    assert list(correlation.columns) == ["col_a", "col_b"]
    assert list(correlation.index) == ["col_a", "col_b"]
    assert correlation.loc["col_a", "col_a"] == pytest.approx(1.0)
    assert correlation.loc["col_b", "col_b"] == pytest.approx(1.0)
    assert correlation.loc["col_a", "col_b"] == pytest.approx(1.0)


def test_build_correlation_matrix__returns_empty_when_no_numeric():
    frame = pd.DataFrame({"letters": ["a", "b", "c"]})

    correlation = build_correlation_matrix(frame, table_name="letters")

    assert correlation.empty
    assert list(correlation.columns) == []
