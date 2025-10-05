from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import iuphar


@pytest.fixture
def sample_input(tmp_path: Path) -> Path:
    path = tmp_path / "output.target_20240101.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL1",
                "GuidetoPHARMACOLOGY": "GTOP1",
                "iuphar_target_id": "T1",
                "iuphar_family_id": "F1",
                "iuphar_type": "Receptor",
                "iuphar_class": "ClassA",
                "iuphar_subclass": "SubClass",
                "iuphar_chain": "ChainA",
                "iuphar_name": "Alpha Beta",
                "gtop_synonyms": "Alpha|Beta|alpha",
                "synonyms": "Alpha|Delta",
                "component_description": '[{"description": "Gamma (extra)"}]',
                "component_synonym_ids": "drop-me",
                "component_type_raw": "unused",
            },
            {
                "target_chembl_id": "CHEMBL2",
                "GuidetoPHARMACOLOGY": "GTOP2",
                "iuphar_target_id": "T2",
                "iuphar_family_id": "F2",
                "iuphar_type": "Enzyme",
                "iuphar_class": "ClassB",
                "iuphar_subclass": "SubClass",
                "iuphar_chain": "ChainB",
                "iuphar_name": "Delta",
                "gtop_synonyms": "",
                "synonyms": "",
                "component_description": "Sigma|Delta",
            },
        ]
    )
    frame.to_csv(path, index=False)
    return path


def test_process_iuphar_targets__transforms_synonyms(tmp_path: Path, sample_input: Path, caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.INFO)

    output_path = iuphar.process_iuphar_targets(sample_input, verbose=True)

    expected_output = sample_input.with_name("IUPHAR.output.target_20240101.csv")
    assert output_path == expected_output
    assert output_path.exists()

    result = pd.read_csv(output_path)
    assert result.columns.tolist() == [
        "target_chembl_id",
        "guidetopharmacology_id",
        "iuphar_target_id",
        "iuphar_family_id",
        "iuphar_type",
        "iuphar_class",
        "iuphar_subclass",
        "iuphar_chain",
        "iuphar_name",
        "iuphar_synonyms",
    ]
    assert result.loc[0, "guidetopharmacology_id"] == "GTOP1"
    assert result.loc[0, "iuphar_synonyms"] == "alpha|beta|delta|gamma|alpha beta"
    assert result.loc[1, "iuphar_synonyms"] == "sigma|delta"

    log_lines = "\n".join(caplog.messages)
    assert "input_rows=2" in log_lines
    assert "output_rows=2" in log_lines
    assert "dropped_columns=2" in log_lines
    assert "synonym_tokens_before=10" in log_lines
    assert "synonym_tokens_after=7" in log_lines


def test_process_iuphar_targets__auto_discovers_latest(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    first = tmp_path / "output.target_20230101.csv"
    second = tmp_path / "output.target_20240202.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL3",
                "GuidetoPHARMACOLOGY": "GTOP3",
                "iuphar_target_id": "T3",
                "iuphar_family_id": "F3",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Sub",
                "iuphar_chain": "Chain",
                "iuphar_name": "Theta",
                "gtop_synonyms": "",
                "synonyms": "",
                "component_description": "",
            }
        ]
    )
    frame.to_csv(first, index=False)
    frame.to_csv(second, index=False)
    monkeypatch.setattr(iuphar, "_DEFAULT_SEARCH_DIR", tmp_path)

    output_path = iuphar.process_iuphar_targets()

    assert output_path == tmp_path / "IUPHAR.output.target_20240202.csv"
    result = pd.read_csv(output_path)
    assert result.loc[0, "guidetopharmacology_id"] == "GTOP3"


def test_process_iuphar_targets__missing_required_columns(tmp_path: Path) -> None:
    path = tmp_path / "output.target_20240505.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL4",
                "GuidetoPHARMACOLOGY": "GTOP4",
                "iuphar_target_id": "T4",
                "iuphar_family_id": "F4",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Sub",
                "iuphar_chain": "Chain",
                "iuphar_name": "Name",
                # Missing gtop_synonyms and component_description
            }
        ]
    )
    frame.to_csv(path, index=False)

    with pytest.raises(iuphar.IUPHARPostProcessingError) as excinfo:
        iuphar.process_iuphar_targets(path)

    message = str(excinfo.value)
    assert "gtop_synonyms" in message
    assert "component_description" in message
