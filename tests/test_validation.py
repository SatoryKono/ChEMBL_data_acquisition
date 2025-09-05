"""Tests for :mod:`library.validation`."""

from __future__ import annotations

import pandas as pd
import pytest

from library import validation


def test_validate_columns_ok() -> None:
    df = pd.DataFrame({"a": [1], "b": [2]})
    validation.validate_columns(df, ["a", "b"])


def test_validate_columns_missing() -> None:
    df = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError):
        validation.validate_columns(df, ["a", "b"])
