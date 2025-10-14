from __future__ import annotations

import json
from datetime import UTC, datetime

import pandas as pd
from scripts import make_assay_postprocessing as cli

from library.postprocessing.assays.schema import ASSAY_SCHEMA


def test_make_assay_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "assays_raw.csv"
    df = pd.DataFrame(
        {
            "assay_chembl_id": ["CHEMBL3", " chembl2"],
            "assay_type": [" primary", "confirmatory"],
            "assay_test_type": ["functional", "BINDING"],
            "description": [" first entry ", "second"],
            "target_chembl_id": [" t1", "t2"],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "assays_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    expected_columns = [
        col for col in ASSAY_SCHEMA.column_order if col in result.columns
    ]
    assert list(result.columns) == expected_columns
    assert result["assay_chembl_id"].tolist() == ["CHEMBL2", "CHEMBL3"]
    assert result["assay_type"].tolist() == ["CONFIRMATORY", "PRIMARY"]
    assert result["assay_test_type"].tolist() == ["BINDING", "FUNCTIONAL"]
    assert result["description"].tolist() == ["second", " first entry "]
    assert result["target_chembl_id"].tolist() == ["T2", "T1"]
    assert result["is_confirmatory"].tolist() == [True, True]
    assert result["pipeline_version"].notna().all()

    report_path = output_path.parent / "assays.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "assays"
    assert payload["metrics"]["output"]["rows"] == len(result)
    extras = payload.get("extras")
    assert extras and extras["output_postprocessed"] == str(output_path)

    date_str = datetime.now(UTC).strftime("%Y%m%d")
    log_path = log_dir / f"make_assays_postprocessing_{date_str}.log"
    assert log_path.exists()

    snapshot = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == snapshot
