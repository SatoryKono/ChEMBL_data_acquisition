from __future__ import annotations

import json
from datetime import UTC, datetime

import pandas as pd

from library.postprocessing.targets.schema import TARGET_SCHEMA
from scripts import make_target_postprocessing as cli


def test_make_target_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "targets_raw.csv"
    df = pd.DataFrame(
        {
            "target_chembl_id": ["t2", " t1"],
            "pref_name": [" Kinase ", "Receptor"],
            "target_type": ["protein", "PROTEIN FAMILY"],
            "protein_classifications": [
                json.dumps(
                    [
                        {
                            "protein_classification": {
                                "pref_name": "Enzyme",
                                "class_level": 1,
                            }
                        },
                        {
                            "protein_classification": {
                                "pref_name": "Kinase family",
                                "class_level": 2,
                            }
                        },
                    ]
                ),
                "",
            ],
            "protein_class_pred_L1": [None, "Other Protein Target"],
            "protein_class_pred_L2": [None, "Miscellaneous"],
            "synonyms": ["alpha;beta", None],
            "chembl_synonyms": [None, "gamma"],
            "gtopdb_synonyms": ["delta", None],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "targets_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    expected_columns = [
        col for col in TARGET_SCHEMA.column_order if col in result.columns
    ]
    assert list(result.columns)[: len(expected_columns)] == expected_columns
    assert {"chembl_synonyms", "gtopdb_synonyms"}.issubset(result.columns)
    assert result["target_chembl_id"].tolist() == ["T1", "T2"]
    assert result["pref_name"].tolist() == ["Receptor", "Kinase"]
    synonyms = result["synonyms"].tolist()
    assert {token.strip() for token in synonyms[0].split(";")} == {"Receptor", "gamma"}
    assert {token.strip() for token in synonyms[1].split(";")} == {"Kinase", "alpha", "beta", "delta"}
    chembl_synonyms = result["chembl_synonyms"].tolist()
    gtopdb_synonyms = result["gtopdb_synonyms"].tolist()
    assert chembl_synonyms[0] == "gamma"
    assert chembl_synonyms[1] == "Kinase"
    assert pd.isna(gtopdb_synonyms[0])
    assert gtopdb_synonyms[1] == "delta"
    assert result["organism"].isna().all()
    assert result["target_class"].tolist() == ["Other Protein Target", "Enzyme"]
    assert result["protein_family"].tolist() == ["Miscellaneous", "Kinase family"]
    assert result["pipeline_version"].notna().all()

    report_path = output_path.parent / "targets.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "targets"
    assert payload["metrics"]["output"]["rows"] == len(result)
    extras = payload.get("extras")
    assert extras and extras["output_postprocessed"] == str(output_path)

    date_str = datetime.now(UTC).strftime("%Y%m%d")
    log_path = log_dir / f"make_targets_postprocessing_{date_str}.log"
    assert log_path.exists()

    snapshot = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == snapshot
