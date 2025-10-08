from __future__ import annotations

from datetime import datetime, timezone
import json

import pandas as pd

from library.postprocess.documents.schema import DOCUMENT_SCHEMA

from scripts import make_document_postprocessing as cli


def test_make_document_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "documents_raw.csv"
    df = pd.DataFrame(
        {
            "document_chembl_id": ["CHEMBL_DOC2", "CHEMBL_DOC1"],
            "title": ["  First Document  ", "Second"],
            "doc_type": ["article", "journal"],
            "year": [2020, 2019],
            "journal": ["  Nature ", "Science"],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "documents_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    expected_columns = [col for col in DOCUMENT_SCHEMA.column_order if col in result.columns]
    assert list(result.columns) == expected_columns
    assert result["document_chembl_id"].tolist() == ["CHEMBL_DOC1", "CHEMBL_DOC2"]
    assert result["title"].tolist() == ["Second", "First Document"]
    assert result["doc_type"].tolist() == ["journal", "article"]
    assert result["publication_year"].tolist() == [2019, 2020]
    assert result["year"].tolist() == [2019, 2020]
    assert result["journal"].tolist() == ["Science", "Nature"]
    assert result["pipeline_version"].notna().all()

    report_path = output_path.parent / "documents.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "documents"
    assert payload["metrics"]["output"]["rows"] == len(result)

    date_str = datetime.now(timezone.utc).strftime("%Y%m%d")
    log_path = log_dir / f"make_document_postprocessing_{date_str}.log"
    assert log_path.exists()

    snapshot = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == snapshot

