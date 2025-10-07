from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.assays import steps


@pytest.mark.unit
def test_normalize_assay_metadata__applies_flags() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["  confirm assay  "],
            "assay_test_type": ["Primary   screen"],
            "assay_format": ["Plate-based "],
        }
    )

    result = steps.normalize_assay_metadata(frame)

    assert result.loc[0, "assay_type"] == "CONFIRM ASSAY"
    assert result.loc[0, "assay_test_type"] == "PRIMARY SCREEN"
    assert result.loc[0, "assay_format"] == "PLATE-BASED"

    custom = steps.normalize_assay_metadata(frame, uppercase_categories=False)
    assert custom.loc[0, "assay_type"] == "confirm assay"
    assert custom.loc[0, "assay_test_type"] == "Primary screen"


@pytest.mark.unit
def test_normalize_assay_metadata__retains_leading_whitespace_when_requested() -> None:
    frame = pd.DataFrame({"assay_type": ["  Primary   screen"], "assay_format": [pd.NA]})

    result = steps.normalize_assay_metadata(
        frame,
        strip_whitespace=False,
        uppercase_categories=False,
    )

    assert result.loc[0, "assay_type"].startswith(" ")
    assert result.loc[0, "assay_type"].endswith("screen")


@pytest.mark.unit
def test_enrich_assay_flags__uses_configured_terms() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": [
                "Primary confirmation",
                "screening",
                pd.NA,
            ]
        }
    )

    result = steps.enrich_assay_flags(
        frame,
        confirmatory_terms=("primary", "confirm"),
        default_flag=False,
    )

    assert result["is_confirmatory"].tolist() == [True, False, False]


@pytest.mark.unit
def test_enrich_assay_flags__missing_column_uses_default_flag() -> None:
    frame = pd.DataFrame({"assay_format": ["plate"]})

    result = steps.enrich_assay_flags(frame, default_flag=True)

    assert result["is_confirmatory"].tolist() == [True]
