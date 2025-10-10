from __future__ import annotations

from codecs import BOM_UTF8
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import iuphar


def test_ensure_required_columns__raises_for_missing_columns() -> None:
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL1",
                "GuidetoPHARMACOLOGY": "GTOP1",
            }
        ]
    )

    with pytest.raises(iuphar.IUPHARPostProcessingError) as excinfo:
        iuphar._ensure_required_columns(frame)

    message = str(excinfo.value)
    assert "gtop_synonyms" in message

    assert "component_description" not in message


def test_clean_brackets__removes_nested_annotations() -> None:
    value = "Alpha (beta) [legacy] gamma"

    assert iuphar._clean_brackets(value) == "Alpha gamma"


def test_collect_synonym_tokens__handles_none_values() -> None:
    row = pd.Series(
        {
            "gtop_synonyms": None,
            "synonyms": None,
            "component_description": None,
            "iuphar_name": "Omega Chain",
        }
    )

    cleaned, deduped = iuphar._collect_synonym_tokens(row)

    assert cleaned == ["omega chain"]
    assert deduped == ["omega chain"]


def test_collect_synonym_tokens__deduplicates_preserving_order() -> None:
    row = pd.Series(
        {
            "gtop_synonyms": "Alpha|Alpha|Beta",
            "synonyms": "beta|Gamma|beta",
            "component_description": "Alpha; Delta",
            "iuphar_name": "Alpha",
        }
    )

    cleaned, deduped = iuphar._collect_synonym_tokens(row)

    assert cleaned.count("alpha") > 1
    assert cleaned.count("beta") > 1
    assert deduped == ["alpha", "beta", "gamma", "delta"]


def test_parse_component_descriptions__supports_json_payloads() -> None:
    payload = '[{"component_description": "Alpha"}, {"description": "Beta"}, {"name": "Gamma"}]'

    result = iuphar._parse_component_descriptions(payload)

    assert result == ["Alpha", "Beta", "Gamma"]


def test_parse_component_descriptions__handles_fallback_delimiters() -> None:
    payload = "Alpha| Beta;Gamma ;; "

    result = iuphar._parse_component_descriptions(payload)

    assert result == ["Alpha", "Beta", "Gamma"]


def test_collect_synonym_tokens__ignores_empty_array_tokens() -> None:
    row = pd.Series(
        {
            "gtop_synonyms": "[]",
            "synonyms": "[]",
            "component_description": "[]",
            "iuphar_name": "Theta",
        }
    )

    cleaned, deduped = iuphar._collect_synonym_tokens(row)

    assert cleaned == ["theta"]
    assert deduped == ["theta"]


def test_prepare_output__orders_columns_and_applies_renames() -> None:
    frame = pd.DataFrame(
        [
            {
                "GuidetoPHARMACOLOGY": "GTOP1",
                "target_chembl_id": "CHEMBL1",
                "iuphar_target_id": "T1",
                "iuphar_family_id": "F1",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Subclass",
                "iuphar_chain": "Chain",
                "iuphar_name": "Name",
            }
        ]
    )

    prepared = iuphar._prepare_output(frame)

    assert list(prepared.columns) == list(iuphar._OUTPUT_COLUMNS)
    assert prepared.loc[0, "guidetopharmacology_id"] == "GTOP1"
    assert prepared.loc[0, "iuphar_synonyms"] == ""


def test_latest_target_file__returns_newest_export(tmp_path: Path) -> None:
    oldest = tmp_path / "output.target_20230101.csv"
    newest = tmp_path / "output.target_20240202.csv"
    for path in (oldest, newest):
        path.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf-8")

    result = iuphar._latest_target_file(tmp_path)

    assert result == newest


def test_process_iuphar_targets__produces_expected_csv(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    input_path = tmp_path / "output.target_20240101.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL123",
                "GuidetoPHARMACOLOGY": "GTOP123",
                "iuphar_target_id": "T-123",
                "iuphar_family_id": "F-10",
                "iuphar_type": "Receptor",
                "iuphar_class": "ClassA",
                "iuphar_subclass": "SubclassA",
                "iuphar_chain": "A",
                "iuphar_name": "Alpha Receptor",
                "gtop_synonyms": "Alpha|Beta|Alpha (v1)",
                "synonyms": "Beta|Gamma|Gamma",
                "component_description": '[{"component_description": "Delta [old]"}, {"name": "ALPHA"}]',
                "component_synonym_ids": "S-1",
                "component_type_raw": "type",
                "component_sequence": "SEQ",
                "component_structures": "STRUCT",
            },
            {
                "target_chembl_id": "CHEMBL999",
                "GuidetoPHARMACOLOGY": "GTOP999",
                "iuphar_target_id": "T-999",
                "iuphar_family_id": "F-99",
                "iuphar_type": "Enzyme",
                "iuphar_class": "ClassB",
                "iuphar_subclass": "SubclassB",
                "iuphar_chain": "B",
                "iuphar_name": "Omega",
                "gtop_synonyms": "[]",
                "synonyms": None,
                "component_description": "[]",
                "component_synonym_ids": None,
                "component_type_raw": None,
                "component_sequence": None,
                "component_structures": None,
            },
        ]
    )
    frame.to_csv(input_path, index=False)

    output_path = iuphar.process_iuphar_targets(input_path)

    assert output_path.name == "IUPHAR.output.target_20240101.csv"
    expected_bytes = (
        snapshot_resource / "iuphar_postprocessing_expected.csv"
    ).read_bytes()
    actual_bytes = output_path.read_bytes()
    # Normalise line endings for cross-platform determinism
    expected_norm = expected_bytes.replace(b"\r\n", b"\n")
    actual_norm = actual_bytes.replace(b"\r\n", b"\n")
    assert actual_norm == expected_norm
    assert not actual_bytes.startswith(BOM_UTF8)
    expected_header = ",".join(iuphar._OUTPUT_COLUMNS) + "\n"
    assert actual_norm.startswith(expected_header.encode("utf-8"))


def test_process_iuphar_targets__normalises_tmp_suffix(tmp_path: Path) -> None:
    input_path = tmp_path / "output.target_20250101.csv.tmp"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL42",
                "GuidetoPHARMACOLOGY": "GTOP42",
                "iuphar_target_id": "T-42",
                "iuphar_family_id": "F-42",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Subclass",
                "iuphar_chain": "A",
                "iuphar_name": "Alpha",
                "gtop_synonyms": "Alpha|Beta",
            }
        ]
    )
    frame.to_csv(input_path, index=False)

    output_path = iuphar.process_iuphar_targets(input_path)

    assert output_path.name == "IUPHAR.output.target_20250101.csv"
    assert output_path.exists()


def test_process_iuphar_targets__fills_missing_component_description(
    tmp_path: Path,
) -> None:
    path = tmp_path / "output.target_20240606.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL5",
                "GuidetoPHARMACOLOGY": "GTOP5",
                "iuphar_target_id": "T5",
                "iuphar_family_id": "F5",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Sub",
                "iuphar_chain": "Chain",
                "iuphar_name": "Name",
                "gtop_synonyms": "Alpha|Beta",
                "synonyms": "Gamma",
                # component_description intentionally omitted
            }
        ]
    )
    frame.to_csv(path, index=False)

    output_path = iuphar.process_iuphar_targets(path)

    assert output_path.exists()
    result = pd.read_csv(output_path)
    assert result.loc[0, "iuphar_synonyms"] == "alpha|beta|gamma|name"


def test_process_iuphar_targets__handles_header_whitespace(tmp_path: Path) -> None:
    path = tmp_path / "output.target_20240707.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL6",
                "GuidetoPHARMACOLOGY": "GTOP6",
                "iuphar_target_id": "T6",
                "iuphar_family_id": "F6",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Sub",
                "iuphar_chain": "Chain",
                "iuphar_name": "Omega",
                "gtop_synonyms": "Omega|Omega (alt)",
                "component_description": "Omega component",
            }
        ]
    )
    frame = frame.rename(columns={col: f"  {col}  " for col in frame.columns})
    frame.to_csv(path, index=False)

    output_path = iuphar.process_iuphar_targets(path)

    assert output_path.exists()
    result = pd.read_csv(output_path)
    assert result.loc[0, "iuphar_synonyms"] == "omega|omega component"


def test_process_iuphar_targets__handles_bom_prefixed_header(tmp_path: Path) -> None:
    path = tmp_path / "output.target_20240808.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL7",
                "GuidetoPHARMACOLOGY": "GTOP7",
                "iuphar_target_id": "T7",
                "iuphar_family_id": "F7",
                "iuphar_type": "Type",
                "iuphar_class": "Class",
                "iuphar_subclass": "Sub",
                "iuphar_chain": "Chain",
                "iuphar_name": "Lambda",
                "gtop_synonyms": "Lambda",
                "component_description": "Lambda",
            }
        ]
    )
    frame.to_csv(path, index=False, encoding="utf-8-sig")

    output_path = iuphar.process_iuphar_targets(path)

    assert output_path.exists()
    result = pd.read_csv(output_path)
    assert result.loc[0, "guidetopharmacology_id"] == "GTOP7"
