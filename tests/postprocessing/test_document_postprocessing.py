from __future__ import annotations

import pandas as pd
import pytest

from library.postprocessing import document


@pytest.mark.unit
def test_load_output_document__coerces_invalid_numeric_ranges(tmp_path):
    columns = list(document.OUTPUT_DTYPE.keys())
    row = {column: "" for column in columns}

    row.update(
        {
            "ChEMBL.document_chembl_id": "DOC0001",
            "ChEMBL.issue": "11-12",
            "ChEMBL.year": "2020",
            "PubMed.PMID": "12345",
        }
    )

    frame = pd.DataFrame([row], columns=columns)
    csv_path = tmp_path / "output.csv"
    frame.to_csv(csv_path, index=False)

    loaded = document._load_output_document(csv_path)

    assert loaded["ChEMBL.issue"].dtype == "Int64"
    assert pd.isna(loaded.at[0, "ChEMBL.issue"])
    assert loaded["ChEMBL.year"].dtype == "Int64"
    assert loaded.at[0, "ChEMBL.year"] == 2020
    assert loaded["PubMed.PMID"].dtype == "Int64"
    assert loaded.at[0, "PubMed.PMID"] == 12345
