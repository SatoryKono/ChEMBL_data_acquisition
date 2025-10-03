"""Tests for :mod:`library.pipelines.common.normalization`."""

from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.common.normalization import (
    normalize_activities,
    normalize_assays,
    normalize_documents,
    normalize_targets,
    normalize_testitems,
)


@pytest.mark.parametrize(
    "func",
    [
        normalize_activities,
        normalize_assays,
        normalize_documents,
        normalize_targets,
        normalize_testitems,
    ],
)
def test_normalization_common(func) -> None:
    """All normalisation helpers perform standard replacements."""
    df = pd.DataFrame(
        {
            "any_id": [" CHEMBL1 "],
            "any_relation": ["<"],
            "any_units": ["5 μM"],
        }
    )
    original = df.copy()
    result = func(df)

    assert result.loc[0, "any_id"] == "CHEMBL1"
    assert result.loc[0, "any_relation"] == "<="
    assert result.loc[0, "any_units"] == "5 uM"

    # Ensure original data is not modified
    assert original.equals(df)
