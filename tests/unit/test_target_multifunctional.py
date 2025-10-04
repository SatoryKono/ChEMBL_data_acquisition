from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing.target.multifunctional import (
    _COLUMNS_TO_REMOVE,
    compute_multifunctional,
)


@pytest.mark.unit
def test_compute_multifunctional__derives_boolean_flag() -> None:
    base_rows = [
        {
            "target_chembl_id": "CHEMBL1",
            "reaction_ec_numbers": "1.2.3.4|1.2.3.5|2.7.11.1",
        },
        {
            "target_chembl_id": "CHEMBL2",
            "reaction_ec_numbers": "3.1.1.1",
        },
        {
            "target_chembl_id": "CHEMBL3",
            "reaction_ec_numbers": "",
        },
    ]
    rows: list[dict[str, str]] = []
    for row in base_rows:
        enriched = row.copy()
        for column in _COLUMNS_TO_REMOVE:
            enriched.setdefault(column, "")
        rows.append(enriched)
    source = pd.DataFrame(rows)
    result = compute_multifunctional(source)
    assert bool(result.loc[0, "multifunctional_enzyme"]) is True
    assert bool(result.loc[1, "multifunctional_enzyme"]) is False
    assert bool(result.loc[2, "multifunctional_enzyme"]) is False
    assert result.loc[0, "reaction_ec_numbers"] == ["1", "2"]
    assert result.loc[2, "reaction_ec_numbers"] == [""]


@pytest.mark.unit
def test_compute_multifunctional__ignores_missing_optional_columns() -> None:
    source = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL1",
                "reaction_ec_numbers": "1.2.3.4|2.7.11.1",
            }
        ]
    )

    result = compute_multifunctional(source)

    assert result.columns.tolist() == [
        "target_chembl_id",
        "reaction_ec_numbers",
        "multifunctional_enzyme",
    ]
    assert result.loc[0, "reaction_ec_numbers"] == ["1", "2"]
    assert bool(result.loc[0, "multifunctional_enzyme"]) is True
