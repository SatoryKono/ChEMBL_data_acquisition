"""Unit tests for :mod:`library.utils.data_correlation`."""

from __future__ import annotations

import pandas as pd
import pytest

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
