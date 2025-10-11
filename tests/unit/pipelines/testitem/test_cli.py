from pathlib import Path

import pandas as pd
import pytest

from library.config import Config
from library.pipelines.testitem import ParentLookupStats, finalize_output


@pytest.mark.unit
def test_finalize_output__working_output_path_normalised(
    tmp_path: Path, cfg: Config
) -> None:
    working_output = tmp_path / ".output.testitem_20240101.csv.tmp"
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding="utf-8")

    rows = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    parent_stats = ParentLookupStats(
        source="skipped",
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )

    exit_code, artifacts = finalize_output(
        [rows],
        cfg=cfg,
        output=working_output,
        parent_stats_supplier=lambda: parent_stats,
        input_csv=input_csv,
        missing_ids=(),
        emit_legacy_artifacts=False,
    )

    expected_dataset = tmp_path / "output.testitem_20240101.csv"
    expected_corr = (
        tmp_path / "output.testitem_20240101_data_correlation_report_table.csv"
    )
    expected_qc = tmp_path / "output.testitem_20240101_quality_report_table.csv"

    assert exit_code == 0
    assert artifacts is not None
    assert artifacts.dataset == expected_dataset
    assert artifacts.correlation_report == expected_corr
    assert artifacts.quality_report == expected_qc
    assert expected_dataset.exists()
    assert expected_corr.exists()
    assert expected_qc.exists()
    assert not working_output.exists()
