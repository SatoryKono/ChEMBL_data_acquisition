"""Tests for :mod:`library.common.pandas_utils`."""

from __future__ import annotations

import pandas as pd

from library.common.pandas_utils import merge_series_prefer_left


def test_merge_series_prefer_left__fills_missing_with_duplicate_index():
    left = pd.Series([None, "kept"], index=["CHEMBL1", "CHEMBL1"])
    right = pd.Series(["fallback"], index=["CHEMBL1"])

    result = merge_series_prefer_left(left, right)

    expected = pd.Series(["fallback", "kept"], index=["CHEMBL1", "CHEMBL1"])
    pd.testing.assert_series_equal(result, expected)


def test_merge_series_prefer_left__preserves_order_for_multiple_matches():
    left = pd.Series([None, None, "base"], index=["A", "A", "B"])
    right = pd.Series(["first", "second", "fallback"], index=["A", "A", "B"])

    result = merge_series_prefer_left(left, right)

    expected = pd.Series(["first", "second", "base"], index=["A", "A", "B"])
    pd.testing.assert_series_equal(result, expected)


def test_merge_series_prefer_left__ignores_extra_right_rows():
    left = pd.Series([None], index=["X"])
    right = pd.Series(["fill", "extra"], index=["X", "X"])

    result = merge_series_prefer_left(left, right)

    expected = pd.Series(["fill"], index=["X"])
    pd.testing.assert_series_equal(result, expected)
