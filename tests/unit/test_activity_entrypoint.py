"""Unit tests for helpers in :mod:`library.cli.entrypoints.activity`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.cli.entrypoints import activity


@pytest.mark.unit
@pytest.mark.parametrize(
    ("dtype", "expected_dtype", "values"),
    [
        ("Float64", "Float64", [1.5, None, 2.0]),
        ("Int64", "Int64", [1, None, 3]),
        ("boolean", "boolean", [True, None, False]),
        ("string", "string", ["foo", None, "bar"]),
    ],
)
def test_coerce_series_dtype__extension_roundtrip(dtype: str, expected_dtype: str, values: list[object]) -> None:
    series = pd.Series(values)

    result = activity._coerce_series_dtype(series, dtype)

    assert str(result.dtype) == expected_dtype
    expected_series = pd.Series(values, dtype=expected_dtype)
    pd.testing.assert_series_equal(result, expected_series)


@pytest.mark.unit
@pytest.mark.parametrize(
    "dtype",
    ["category", "invalid-type"],
)
def test_coerce_series_dtype__fallback_to_string(dtype: str) -> None:
    series = pd.Series(["1", "2", None])

    result = activity._coerce_series_dtype(series, dtype)

    assert result.dtype == object
    assert result.tolist() == ["1", "2", "None"]


@pytest.mark.unit
@pytest.mark.parametrize(
    ("value", "default", "expected"),
    [
        (5, 0, 5),
        ("7", 0, 7),
        (None, 3, 3),
        ("not-a-number", 2, 2),
        (["1"], 4, 4),
    ],
)
def test_safe_int__conversion_cases(value: object, default: int, expected: int) -> None:
    assert activity._safe_int(value, default) == expected
