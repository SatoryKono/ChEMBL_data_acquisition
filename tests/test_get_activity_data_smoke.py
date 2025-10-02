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

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        int_ids = [int(i) for i in ids]
        return pd.DataFrame(
            {
                "activity_id": int_ids,
                "molecule_chembl_id": ["CHEMBL1" for _ in int_ids],
                "assay_chembl_id": ["CHEMBL_ASSAY" for _ in int_ids],
                "target_id": ["CHEMBL2" for _ in int_ids],
                "standard_type": ["IC50" for _ in int_ids],
                "standard_value": [1.0 for _ in int_ids],
                "pA_value": [5.0 for _ in int_ids],
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)
    output_csv = tmp_path / "out.csv"
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text("activity:\n  column: activity_id\n")
    monkeypatch.setattr(get_activity_data, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(get_activity_data, "write_meta_yaml", lambda **__: None)
    exit_code = get_activity_data.main(
        [
            "--input",
            str(input_csv),
            "--final-out",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(cfg_file),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty


def test_get_activity_data_streaming_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "activity_ids.csv"
    ids = pd.DataFrame({"activity_id": [str(i) for i in range(1, 10)]})
    ids.to_csv(input_csv, index=False)

    calls: list[list[str]] = []

    def fake_get(ids, cfg, client, chunk_size, timeout, **kwargs):
        calls.append(list(ids))
        int_ids = [int(i) for i in ids]
        return pd.DataFrame(
            {
                "activity_id": int_ids,
                "molecule_chembl_id": ["CHEMBL1" for _ in int_ids],
                "assay_chembl_id": ["CHEMBL2" for _ in int_ids],
                "standard_type": ["IC50" for _ in int_ids],
                "standard_value": [1.0 for _ in int_ids],
            }
        )

    monkeypatch.setattr(cl, "get_activities", fake_get)
    output_csv = tmp_path / "stream.csv"
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text(
        "\n".join(
            [
                "activity:",
                "  column: activity_id",
                "  batch_size: 3",
            ]
        )
    )
    monkeypatch.setattr(get_activity_data, "file_sha256", lambda p: "cafebabe")
    monkeypatch.setattr(get_activity_data, "write_meta_yaml", lambda **__: None)
    exit_code = get_activity_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(cfg_file),
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert len(df) == len(ids)
    assert len(calls) >= 3
    for call in calls:
        assert len(call) <= 3
