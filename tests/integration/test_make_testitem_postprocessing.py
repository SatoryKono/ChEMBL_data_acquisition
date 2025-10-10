from __future__ import annotations

import json
from datetime import UTC, datetime

import pandas as pd

from library.postprocessing.testitem.schema import TESTITEM_SCHEMA
from scripts import make_testitem_postprocessing as cli


def test_make_testitem_postprocessing__end_to_end(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs"
    monkeypatch.setenv("CHEMBL_POSTPROCESS_LOG_DIR", str(log_dir))

    input_path = tmp_path / "testitems_raw.csv"
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["chembl2", " CHEMBL1"],
            "parent_molecule_chembl_id": [" chembl3", None],
            "pref_name": [" Test Item ", "Compound"],
            "natural_product": ["YES", ""],
            "prodrug": [0, 1],
            "polymer_flag": [None, "true"],
            "pubchem_cid": ["123", ""],
            "timestamp_utc": ["2021-01-01T00:00:00Z", None],
        }
    )
    df.to_csv(input_path, index=False)

    output_path = tmp_path / "testitems_processed.csv"
    exit_code = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert exit_code == 0

    result = pd.read_csv(output_path)
    assert list(result.columns) == list(TESTITEM_SCHEMA.column_order)
    assert result["molecule_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    parents = result["parent_molecule_chembl_id"].tolist()
    assert pd.isna(parents[0])
    assert parents[1] == "CHEMBL3"
    assert result["pref_name"].tolist() == ["Compound", "Test Item"]
    natural = result["natural_product"].astype("boolean")
    polymer = result["polymer_flag"].astype("boolean")
    assert pd.isna(natural.iloc[0])
    assert bool(natural.iloc[1]) is True
    assert bool(polymer.iloc[0]) is True
    assert pd.isna(polymer.iloc[1])
    assert result["prodrug"].astype("boolean").tolist() == [True, False]
    assert result["pipeline_version"].notna().all()

    report_path = output_path.parent / "testitems.postprocess.report.json"
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["table"] == "testitems"
    assert payload["metrics"]["output"]["rows"] == len(result)

    date_str = datetime.now(UTC).strftime("%Y%m%d")
    log_path = log_dir / f"make_testitem_postprocessing_{date_str}.log"
    assert log_path.exists()

    snapshot = output_path.read_bytes()
    second_exit = cli.main(["--input", str(input_path), "--output", str(output_path)])
    assert second_exit == 0
    assert output_path.read_bytes() == snapshot
