from __future__ import annotations

from datetime import datetime, timezone
import json

import pandas as pd

from library.postprocessing.activities.schema import ACTIVITY_SCHEMA

from scripts import make_activity_postprocessing as cli


def test_make_activity_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "activities_raw.csv"
    df = pd.DataFrame(
        {
            "activity_id": [2, 1],
            "molecule_chembl_id": ["CHEMBL2", "CHEMBL1"],
            "assay_chembl_id": ["ASSAY2", "ASSAY1"],
            "standard_type": ["IC50", "Ki"],
            "standard_relation": ["=", "<"],
            "standard_value": [10.0, 5.0],
            "standard_units": ["nM", "nM"],
            "data_validity_comment": ["Valid result", None],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "activities_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    expected_columns = [col for col in ACTIVITY_SCHEMA.column_order if col in result.columns]
    assert list(result.columns) == expected_columns
    assert result["molecule_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    assert result["standard_units"].tolist() == ["NM", "NM"]
    assert result["quality_flag"].tolist() == [False, True]
    assert result["pipeline_version"].notna().all()

    # Ensure deterministic ordering and schema compliance
    expected_sort = result.sort_values(
        ["molecule_chembl_id", "assay_chembl_id", "activity_id"],
        kind="mergesort",
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(result, expected_sort)

    # Metrics report is generated next to the output
    report_path = output_path.parent / "activities.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "activities"
    assert payload["output_path"] == str(output_path)
    assert payload["metrics"]["output"]["rows"] == len(result)

    # Log file is created with the expected naming convention
    date_str = datetime.now(timezone.utc).strftime("%Y%m%d")
    log_path = log_dir / f"make_activity_postprocessing_{date_str}.log"
    assert log_path.exists()

    # Idempotent reruns keep the output stable
    initial_digest = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == initial_digest

