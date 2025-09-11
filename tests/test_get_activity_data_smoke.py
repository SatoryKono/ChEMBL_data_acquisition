from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library import chembl_library as cl
from scripts import get_activity_data


def test_get_activity_data_smoke(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = Path("tests/data/activity_ids_small.csv")

    def fake_get(ids, cfg, chunk_size, timeout):
        int_ids = [int(i) for i in ids]
        return pd.DataFrame(
            {
                "activity_id": int_ids,
                "testitem_id": ["CHEMBL1" for _ in int_ids],
                "target_id": ["CHEMBL2" for _ in int_ids],
                "standard_type": ["IC50" for _ in int_ids],
                "standard_value": [1.0 for _ in int_ids],
                "pA_value": [5.0 for _ in int_ids],
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)
    output_csv = tmp_path / "out.csv"
    exit_code = get_activity_data.main(
        ["--input", str(input_csv), "--output", str(output_csv), "--log-level", "ERROR"]
    )
    assert exit_code == 0
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty
