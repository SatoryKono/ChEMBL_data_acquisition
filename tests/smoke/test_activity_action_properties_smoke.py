from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from library.config import Config
from tests.pipelines.test_activity_extraction import DummyChemblClient  # type: ignore


def test_activity_action_properties_cli(tmp_path: Path, monkeypatch) -> None:
    from scripts import get_activity_data

    input_csv = tmp_path / "activities.csv"
    input_csv.write_text("activity_id\nA1\nA2\n", encoding="utf-8")
    output_csv = tmp_path / "result.csv"

    df = pd.DataFrame(
        [
            {
                "activity_id": "A1",
                "molecule_chembl_id": "CHEMBL1",
                "assay_chembl_id": "AS1",
                "standard_value": 1.0,
                "standard_type": "IC50",
                "activity_comment": "positive allosteric modulator",
                "assay_description": "Calcium flux",
            },
            {
                "activity_id": "A2",
                "molecule_chembl_id": "CHEMBL2",
                "assay_chembl_id": "AS1",
                "standard_value": 2.0,
                "standard_type": "IC50",
                "activity_comment": "negative allosteric modulator",
                "assay_description": "Binding displacement",
            },
        ]
    )

    monkeypatch.setattr("scripts.get_activity_data.ChemblClient", DummyChemblClient)
    monkeypatch.setattr(
        "scripts.get_activity_data.cl.get_activities", lambda *_, **__: df
    )
    quality_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def _capture_quality(*args: object, **kwargs: object) -> None:
        quality_calls.append((args, kwargs))

    monkeypatch.setattr(
        "scripts.get_activity_data.analyze_table_quality", _capture_quality
    )
    monkeypatch.setattr("scripts.get_activity_data.write_meta_yaml", lambda **__: None)
    monkeypatch.setattr("scripts.get_activity_data.file_sha256", lambda _: "deadbeef")

    cfg = Config()
    cfg.io.output_dir = tmp_path / "output"
    cfg.io.cache_dir = tmp_path / "cache"
    cfg.activity.column = "activity_id"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(cfg.model_dump(mode="json"), sort_keys=False),
        encoding="utf-8",
    )

    exit_code = get_activity_data.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
        ]
    )
    assert exit_code == 0

    assert quality_calls, "quality profiler should be invoked"
    quality_args, quality_kwargs = quality_calls[-1]
    assert Path(quality_args[0]) == output_csv
    assert Path(quality_kwargs.get("destination_dir")).resolve() == output_csv.parent.resolve()

    result = pd.read_csv(output_csv)
    assert set(["action_type", "activity_properties", "properties_hash"]).issubset(
        result.columns
    )
    assert result["action_type"].tolist() == ["inhibition", "inhibition"]

    payloads = result["activity_properties"].map(json.loads).tolist()
    assert payloads[0]["effect_features"]["positive"] is True
    assert payloads[1]["effect_features"]["negative"] is True
    assert all(
        isinstance(value, str) and len(value) == 64
        for value in result["properties_hash"].tolist()
    )
