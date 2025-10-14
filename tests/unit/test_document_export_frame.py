import pandas as pd
import pytest
from scripts import get_document_data


@pytest.mark.unit
def test_prepare_export_frame__tuple_export_columns(monkeypatch):
    original_columns = list(get_document_data._EXPORT_COLUMNS)
    dataframe = pd.DataFrame({"ChEMBL.doi": ["10.1000/example"]})

    monkeypatch.setattr(
        get_document_data,
        "_EXPORT_COLUMNS",
        tuple(original_columns),
        raising=False,
    )

    export_frame = get_document_data._prepare_export_frame(dataframe)

    assert list(export_frame.columns) == original_columns
    assert export_frame.loc[0, "ChEMBL.doi"] == "10.1000/example"
    assert export_frame.loc[0, "PubMed.PMID"] == ""
