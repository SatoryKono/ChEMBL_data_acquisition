"""Tests for :mod:`library.common.pandas_utils`."""

from __future__ import annotations

from collections import defaultdict

import pandas as pd
import pytest

pytest.importorskip("hypothesis")

from hypothesis import given, settings, strategies as st

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


def _manual_merge_series(left: pd.Series, right: pd.Series) -> pd.Series:
    if left.empty:
        return left.copy()

    right_map: dict[object, list[object]] = defaultdict(list)
    for label, value in zip(right.index, right.tolist()):
        right_map[label].append(value)

    counters: defaultdict[object, int] = defaultdict(int)
    result_values: list[object] = []
    for label, value in zip(left.index, left.tolist()):
        position = counters[label]
        counters[label] += 1
        if pd.isna(value):
            replacements = right_map.get(label, [])
            if position < len(replacements):
                result_values.append(replacements[position])
                continue
        result_values.append(value)

    return pd.Series(result_values, index=left.index, dtype=object)


_VALUE_STRATEGY = st.one_of(
    st.none(),
    st.just(pd.NA),
    st.integers(-5, 5),
    st.text(alphabet=list("abcXYZ"), min_size=0, max_size=5),
)
_INDEX_STRATEGY = st.integers(-2, 3)
_PAIR_LIST_STRATEGY = st.lists(
    st.tuples(_VALUE_STRATEGY, _INDEX_STRATEGY),
    min_size=0,
    max_size=6,
)


@settings(max_examples=50, derandomize=True, deadline=None)
@given(left_pairs=_PAIR_LIST_STRATEGY, right_pairs=_PAIR_LIST_STRATEGY)
def test_merge_series_prefer_left__matches_manual_alignment(
    left_pairs: list[tuple[object, object]],
    right_pairs: list[tuple[object, object]],
) -> None:
    left_values = [value for value, _ in left_pairs]
    left_index = [label for _, label in left_pairs]
    right_values = [value for value, _ in right_pairs]
    right_index = [label for _, label in right_pairs]

    left_series = pd.Series(left_values, index=left_index or None, dtype=object)
    right_series = pd.Series(right_values, index=right_index or None, dtype=object)

    result = merge_series_prefer_left(left_series, right_series)
    expected = _manual_merge_series(left_series, right_series)

    pd.testing.assert_series_equal(result, expected, check_dtype=False)
