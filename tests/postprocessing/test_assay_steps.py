from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.assays.steps import normalize_assay_metadata


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
