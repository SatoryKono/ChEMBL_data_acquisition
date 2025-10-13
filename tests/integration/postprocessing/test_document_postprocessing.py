from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from library.config import IoCfg
from library.pipelines.document import postprocessing as stage_postprocessing
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


@pytest.mark.integration
def test_postprocess_file__forwards_custom_reference(tmp_path, monkeypatch):
    input_path = tmp_path / "documents_raw.csv"
    pd.DataFrame({"ChEMBL.document_chembl_id": ["DOC0001"]}).to_csv(
        input_path, index=False
    )

    captured: dict[str, object] = {}

    def _stub_postprocess_documents(df, **kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            {
                "completed": pd.Series([False], dtype="boolean"),
                "document_chembl_id": ["DOC0001"],
            }
        )

    monkeypatch.setattr(
        stage_postprocessing, "postprocess_documents", _stub_postprocess_documents
    )

    output_path = tmp_path / "processed.csv"
    reference_override = tmp_path / "reference.csv"
    resources = SimpleNamespace(dictionary_dir=tmp_path / "resources")

    cfg = IoCfg()
    destination = stage_postprocessing.postprocess_file(
        input_path,
        output_path,
        cfg=cfg,
        ref_document_path=reference_override,
        resources=resources,
    )

    assert captured["ref_document_path"] == reference_override
    assert captured["resources"] is resources

    assert destination == output_path
    assert destination.exists()
    result = pd.read_csv(destination)
    assert list(result.columns) == ["completed", "document_chembl_id"]
    assert result.loc[0, "document_chembl_id"] == "DOC0001"
