"""Tests for :mod:`library.validation`."""

from __future__ import annotations

import pandas as pd
import pytest

from library import validation
from schemas import normalize_testitems


def test_validate_columns_ok() -> None:
    df = pd.DataFrame({"a": [1], "b": [2]})
    validation.validate_columns(df, ["a", "b"])


def test_validate_columns_missing() -> None:
    df = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError):
        validation.validate_columns(df, ["a", "b"])


def _activities_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "activity_id": ["1", "2"],
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "assay_chembl_id": ["ASSAY1", "ASSAY2"],
            "standard_value": [1.0, -1.0],
            "standard_type": ["IC50", "IC50"],
        }
    )


def test_validate_activities_returns_dataframe() -> None:
    df = _activities_frame()

    result = validation.validate_activities(df)

    assert isinstance(result, pd.DataFrame)
    assert result.index.tolist() == [0]
    assert result["standard_value"].tolist() == [1.0]
    assert df.index.tolist() == [0, 1]


def test_validate_activities_result_object() -> None:
    df = _activities_frame()

    result = validation.validate_activities(df, return_result=True)

    assert isinstance(result, validation.ValidationResult)
    assert result.has_failures is True
    assert result.failure_cases["index"].dropna().tolist() == [1]
    assert result.data.index.tolist() == [0]


def test_validate_assays_filters_invalid_rows() -> None:
    df = pd.DataFrame(
        {
            "assay_chembl_id": ["A1", 2],
            "document_chembl_id": ["D1", 3],
            "target_chembl_id": ["T1", "T2"],
        }
    )

    filtered = validation.validate_assays(df)

    assert isinstance(filtered, pd.DataFrame)
    assert filtered.index.tolist() == [0]

    result = validation.validate_assays(df, return_result=True)

    assert isinstance(result, validation.ValidationResult)
    assert result.has_failures is True
    assert result.data.index.tolist() == [0]
    assert sorted(set(result.failure_cases["index"].dropna())) == [1]


def test_validate_testitems_filters_invalid_rows() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", 2],
            "pref_name": ["Name", "Other"],
        }
    )

    filtered = validation.validate_testitems(df)

    assert isinstance(filtered, pd.DataFrame)
    assert filtered.index.tolist() == [0]

    result = validation.validate_testitems(df, return_result=True)

    assert isinstance(result, validation.ValidationResult)
    assert result.has_failures is True
    assert result.data.index.tolist() == [0]
    assert sorted(set(result.failure_cases["index"].dropna())) == [1]


def test_validate_testitems_all_valid() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "pref_name": ["Name"],
        }
    )

    result = validation.validate_testitems(df, return_result=True)

    assert isinstance(result, validation.ValidationResult)
    assert result.has_failures is False
    pd.testing.assert_frame_equal(result.data.reset_index(drop=True), df)


def test_validate_testitems_numeric_ids_preserved() -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": [123456],
            "pref_name": ["Name"],
        }
    )

    normalized = normalize_testitems(df)

    assert normalized["molecule_chembl_id"].dtype == pd.StringDtype()

    result = validation.validate_testitems(normalized, return_result=True)

    assert isinstance(result, validation.ValidationResult)
    assert result.has_failures is False
    assert result.data["molecule_chembl_id"].tolist() == ["123456"]
