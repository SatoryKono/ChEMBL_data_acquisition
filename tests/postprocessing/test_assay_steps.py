from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.pipeline.assays.steps import (
    enrich_assay_flags,
    finalize_assay_records,
    normalize_assay_metadata,
)


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_normalize_assay_metadata__default_normalization() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["  primary  screen  "],
            "assay_test_type": [" confirmatory   "],
            "assay_format": [" cell based"],
        }
    )

    result = normalize_assay_metadata(frame)

    assert result.loc[0, "assay_type"] == "PRIMARY SCREEN"
    assert result.loc[0, "assay_test_type"] == "CONFIRMATORY"
    assert result.loc[0, "assay_format"] == "CELL BASED"
    assert str(result["assay_type"].dtype) == "string"


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_normalize_assay_metadata__respects_disabled_flags() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["  primary  screen  "],
            "assay_test_type": [" confirmatory   "],
            "assay_format": [" cell based"],
        }
    )

    result = normalize_assay_metadata(
        frame,
        uppercase_categories=False,
        strip_whitespace=False,
        collapse_internal_whitespace=False,
    )

    assert result.loc[0, "assay_type"] == "  primary  screen  "
    assert result.loc[0, "assay_test_type"] == " confirmatory   "
    assert result.loc[0, "assay_format"] == " cell based"
    assert str(result["assay_type"].dtype) == "string"


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_normalize_assay_metadata__ignores_unknown_parameters() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["primary"],
            "assay_format": ["cell"],
            "assay_test_type": ["screen"],
        }
    )

    result = normalize_assay_metadata(
        frame, uppercase_categories=False, unsupported_flag=True
    )

    expected = frame.copy()
    for column in expected.columns:
        expected[column] = expected[column].astype("string")

    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_enrich_assay_flags__applies_custom_terms() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["Primary screen", "Counter screen"],
        }
    )

    result = enrich_assay_flags(
        frame,
        confirmatory_terms=["primary", "confirm"],
        default_flag=False,
    )

    assert result["is_confirmatory"].tolist() == [True, False]


@pytest.mark.postprocessing
@pytest.mark.usefixtures("deterministic_env")
def test_finalize_assay_records__normalizes_identifiers() -> None:
    frame = pd.DataFrame(
        {
            "assay_chembl_id": ["  chEMBL123  "],
            "assay_type": ["binding"],
            "assay_test_type": ["confirmatory"],
            "description": ["Example"],
            "target_chembl_id": [" chembl999 "],
        }
    )

    result = finalize_assay_records(frame, normalize_identifiers=True)

    assert result.loc[0, "assay_chembl_id"] == "CHEMBL123"
    assert result.loc[0, "target_chembl_id"] == "CHEMBL999"
