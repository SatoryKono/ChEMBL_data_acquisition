
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import target


@pytest.mark.unit
def test_split_pipes__handles_empty_and_whitespace():
    assert target._split_pipes("") == []
    assert target._split_pipes("a|b") == ["a", "b"]
    assert target._split_pipes(" a | | b ") == ["a", "b"]


@pytest.mark.unit
def test_make_triples__pads_shorter_lists():
    triples = target._make_triples(["alpha", "beta"], ["ID1"], ["syn1", "syn2", "syn3"])
    assert triples == [
        {"name": "alpha", "id": "ID1", "synonym": "syn1"},
        {"name": "beta", "id": None, "synonym": "syn2"},
        {"name": None, "id": None, "synonym": "syn3"},
    ]


@pytest.mark.unit
def test_syn_expand__produces_variants():
    assert target._syn_expand("PDE3A") == ["pde3a", "3a"]
    assert target._syn_expand("PLD2A") == ["pld2a", "2a"]


@pytest.mark.unit
def test_tokenize_synonym__pde3a_alpha():
    assert target._tokenize_synonym("PDE3A:alpha") == ["pde3a", "3a", "alpha"]


def _build_sample_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "isoform_synonyms": "PDE3A:Alpha|Beta Variant",
                "isoform_names": "Alpha|Beta",
                "isoform_ids": "EnSp1|ENSP2",
                "uniprot_id_primary": "U1",
                "target_chembl_id": "CHEMBL1",
            }
        ]
    )


@pytest.mark.unit
def test_transform_frame__normalises_names_preserves_ids():
    frame = _build_sample_frame()
    result, _ = target._transform_frame(frame)

    assert all(name == name.lower() for name in result["name"])
    assert "EnSp1" in result["id"].tolist()


@pytest.mark.unit
def test_transform_frame__filters_placeholder_names():
    frame = pd.DataFrame(
        [
            {
                "isoform_synonyms": "SynA|None|n/a",
                "isoform_names": "Alpha|n/a|none|",
                "isoform_ids": "ID1|ID2|ID3",
                "uniprot_id_primary": "U1",
                "target_chembl_id": "CHEMBL1",
            }
        ]
    )
    result, _ = target._transform_frame(frame)

    assert "n/a" not in result["name"].values
    assert "none" not in result["name"].values
    assert "" not in result["name"].values


@pytest.mark.unit
def test_dedup_stage1__matches_table_distinct():
    combined = pd.DataFrame(
        [
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U2",
            },
        ]
    )
    expected = combined.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"], keep="first"
    )
    assert target._dedupe_stage1(combined).equals(expected)


@pytest.mark.unit
def test_stable_sort__preserves_relative_order():
    df = pd.DataFrame(
        [
            {"uniprot_id_primary": "U2", "id": "B", "marker": 1},
            {"uniprot_id_primary": "U1", "id": "B", "marker": 2},
            {"uniprot_id_primary": "U1", "id": "A", "marker": 3},
            {"uniprot_id_primary": "U1", "id": "A", "marker": 4},
        ]
    )
    sorted_df = target._stable_sort(df)
    assert sorted_df[["uniprot_id_primary", "id", "marker"]].to_dict("records") == [
        {"uniprot_id_primary": "U1", "id": "A", "marker": 3},
        {"uniprot_id_primary": "U1", "id": "A", "marker": 4},
        {"uniprot_id_primary": "U1", "id": "B", "marker": 2},
        {"uniprot_id_primary": "U2", "id": "B", "marker": 1},
    ]


@pytest.mark.unit
def test_dedup_stage2__respects_sorted_order():
    df = pd.DataFrame(
        [
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U2",
            },
            {
                "id": "ID1",
                "name": "beta",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
        ]
    )
    deduped = target._dedupe_stage2(df)
    assert deduped.iloc[0]["uniprot_id_primary"] == "U1"
    assert len(deduped) == 2


@pytest.mark.unit
def test_dedup_stage3__drops_duplicates_by_id_name():
    df = pd.DataFrame(
        [
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
            {
                "id": "ID1",
                "name": "alpha",
                "target_chembl_id": "CHEMBL2",
                "uniprot_id_primary": "U2",
            },
            {
                "id": "ID2",
                "name": "beta",
                "target_chembl_id": "CHEMBL1",
                "uniprot_id_primary": "U1",
            },
        ]
    )
    deduped = target._dedupe_stage3(df)
    assert len(deduped) == 2
    assert deduped.iloc[0]["target_chembl_id"] == "CHEMBL1"


@pytest.mark.unit
def test_transform_frame__column_order_exact():
    frame = _build_sample_frame()
    result, _ = target._transform_frame(frame)
    assert list(result.columns) == target.OUTPUT_COLUMNS


@pytest.mark.unit
def test_process_targets__writes_with_isoform_prefix(tmp_path: Path):
    input_path = tmp_path / "output.target_20250101.csv"
    frame = _build_sample_frame()
    frame.to_csv(input_path, index=False, encoding="utf-8")

    output_path = target.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.output.target_20250101.csv"
    assert output_path.parent == tmp_path
    written = pd.read_csv(output_path)
    assert list(written.columns) == target.OUTPUT_COLUMNS


@pytest.mark.unit
def test_process_targets__auto_discovers_latest_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    first = tmp_path / "output.target_20240101.csv"
    second = tmp_path / "output.target_20240201.csv"
    frame = _build_sample_frame()
    frame.to_csv(first, index=False, encoding="utf-8")
    frame.to_csv(second, index=False, encoding="utf-8")
    os.utime(first, (first.stat().st_atime, first.stat().st_mtime - 100))

    monkeypatch.setattr(target, "DEFAULT_INPUT_DIR", tmp_path)
    output_path = target.process_targets(output_csv=tmp_path / "result.csv", verbose=False)

    assert output_path == tmp_path / "result.csv"
    written = pd.read_csv(output_path)
    assert not written.empty


@pytest.mark.unit
def test_transform_frame__handles_mismatched_lengths():
    frame = pd.DataFrame(
        [
            {
                "isoform_synonyms": "SynA|SynB|SynC",
                "isoform_names": "Alpha|Beta",
                "isoform_ids": "ID1",
                "uniprot_id_primary": "U1",
                "target_chembl_id": "CHEMBL1",
            }
        ]
    )
    result, _ = target._transform_frame(frame)

    beta_rows = result[result["name"] == "beta"]
    assert not beta_rows.empty
    assert beta_rows["id"].isna().all()
