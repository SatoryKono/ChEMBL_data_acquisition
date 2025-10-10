from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from library.postprocessing.target.cellularity import add_cellularity_smart


@pytest.mark.unit
def test_add_cellularity_smart__skips_fetch_for_missing_tax_ids() -> None:
    def _fetch_stub(tax_id: object, email: str | None) -> list[str]:
        raise AssertionError("fetcher should not be invoked for missing tax IDs")

    frame = pd.DataFrame(
        {
            "tax_id": pd.Series([pd.NA, np.nan], dtype="object"),
            "superkingdom": pd.Series([None, None], dtype="object"),
            "phylum": pd.Series([None, None], dtype="object"),
        }
    )

    result = add_cellularity_smart(
        frame,
        tax_id_column="tax_id",
        superkingdom_column="superkingdom",
        phylum_column="phylum",
        fetcher=_fetch_stub,
    )

    assert result["cellularity"].tolist() == ["ambiguous", "ambiguous"]
