from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.assay.postprocessing import postprocess_assays


@pytest.mark.unit
def test_postprocess_assays__empty_without_columns() -> None:
    """Ensure post-processing handles a completely column-less frame."""

    df = pd.DataFrame()

    result = postprocess_assays(df)

    assert result.empty
    assert set(result.columns) == {
        "document_chembl_id",
        "target_chembl_id",
        "assay_with_same_target",
    }
    assert str(result["assay_with_same_target"].dtype) == "Int64"


@pytest.mark.unit
def test_postprocess_assays__empty_with_required_columns() -> None:
    """Ensure post-processing preserves required columns on empty input."""

    df = pd.DataFrame(columns=["document_chembl_id", "target_chembl_id"])

    result = postprocess_assays(df)

    assert result.empty
    assert list(result.columns) == [
        "document_chembl_id",
        "target_chembl_id",
        "assay_with_same_target",
    ]
    assert str(result["assay_with_same_target"].dtype) == "Int64"
