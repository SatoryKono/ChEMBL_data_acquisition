from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from scripts import get_target_data

from library.postprocessing.target.steps import TargetData


def test_get_target_data_cli(tmp_path, monkeypatch) -> None:
    dataset = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "pref_name": ["Target 1", "Target 2"],
            "target_type": ["protein", "protein"],
            "organism": ["Human", "Human"],
            "uniprot_id": ["P12345", "P67890"],
            "gene_symbol": ["GENE1", "GENE2"],
            "target_class": ["ClassA", "ClassB"],
            "protein_family": ["FamilyA", "FamilyB"],
            "synonyms": ["Syn1", "Syn2"],
        }
    ).astype("string")

    quality_report = pd.DataFrame({"column": ["target_chembl_id"], "non_null": [1]}).astype("string")
    correlation_report = pd.DataFrame()
    qc_summary = {"row_count": 2, "column_count": 9, "non_null_ratio": 1.0}
    target_data = TargetData(
        dataset=dataset,
        quality_report=quality_report,
        correlation_report=correlation_report,
        qc_summary=qc_summary,
    )

    monkeypatch.setattr(
        "library.postprocessing.target.steps.fetch_normalize_target",
        lambda limit: dataset,
    )
    monkeypatch.setattr(
        "library.postprocessing.target.steps.generate_target_reports",
        lambda frame: target_data,
    )

    exit_code = get_target_data.main(
        [
            "all",
            "--limit",
            "10",
            "--date-tag",
            "20250101",
            "--output-dir",
            str(tmp_path),
            "--log-level",
            "debug",
            "--input-dir",
            str(tmp_path / "input"),
        ]
    )

    assert exit_code == 0

    expected_stem = tmp_path / "output.target_20250101"
    dataset_path = Path(f"{expected_stem}.csv")
    quality_path = Path(f"{expected_stem}_quality_report_table.csv")
    correlation_path = Path(f"{expected_stem}_data_correlation_report_table.csv")
    metadata_path = Path(f"{expected_stem}.meta.yaml")

    assert dataset_path.exists()
    assert quality_path.exists()
    assert correlation_path.exists()
    assert metadata_path.exists()


def test_parse_args_accepts_logging_and_input_dir(tmp_path) -> None:
    namespace, parsed = get_target_data.parse_args(
        [
            "--limit",
            "5",
            "--log-level",
            "debug",
            "--input-dir",
            str(tmp_path),
        ]
    )

    assert parsed.log_level == "DEBUG"
    assert parsed.input_dir == tmp_path
    assert namespace.log_level == "DEBUG"
    assert namespace.input_dir == tmp_path


def test_parse_args_rejects_unknown_log_level() -> None:
    with pytest.raises(SystemExit):
        get_target_data.parse_args(["--log-level", "verbose"])
