from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import target


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("", []),
        ("a|b", ["a", "b"]),
        (" a | | b ", ["a", "b"]),
    ],
)
def test_split_pipes__cases(value: str, expected: list[str]) -> None:
    assert target._split_pipes(value) == expected


def test_make_triples__pads_missing_values() -> None:
    names = ["alpha", "beta"]
    ids = ["ID1"]
    synonyms = ["syn1", "syn2", "syn3"]
    triples = target._make_triples(names, ids, synonyms)
    assert triples == [
        {"name": "alpha", "id": "ID1", "synonym": "syn1"},
        {"name": "beta", "id": None, "synonym": "syn2"},
        {"name": None, "id": None, "synonym": "syn3"},
    ]


def test_syn_expand__removes_pde_and_pld() -> None:
    assert target._syn_expand("PDE3A") == ["pde3a", "3a"]
    assert target._syn_expand("PLDALPHA") == ["pldalpha", "alpha"]


def test_tokenize_synonym__pde3a_alpha() -> None:
    tokens = target._tokenize_synonym("PDE3A:alpha")
    assert tokens == ["pde3a", "3a", "alpha"]


@pytest.fixture()
def minimal_isoform_frame(snapshot_resource: Path) -> pd.DataFrame:
    source = snapshot_resource / "targets_isoform_minimal.csv"
    return pd.read_csv(source, dtype=str, keep_default_na=False, na_filter=False)


def test_lowercase_applied_only_to_names_and_synonyms(
    minimal_isoform_frame: pd.DataFrame,
) -> None:
    processed, _ = target._postprocess_frame(minimal_isoform_frame)
    ids = [value for value in processed["id"].tolist() if value]
    assert "IsoA" in ids
    assert "isoa" not in ids
    assert all(name == name.lower() for name in processed["name"])


def test_invalid_names_filtered(minimal_isoform_frame: pd.DataFrame) -> None:
    processed, _ = target._postprocess_frame(minimal_isoform_frame)
    assert "n/a" not in processed["name"].values
    assert "none" not in processed["name"].values


def test_union_and_dedup_stage1_equivalence() -> None:
    df = pd.DataFrame(
        {
            "id": ["A", "A", "B"],
            "uniprot_id_primary": ["P1", "P1", "P2"],
            "target_chembl_id": ["C1", "C1", "C2"],
            "name": ["alpha", "alpha", "beta"],
        }
    )
    expected = df.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"], keep="first"
    ).reset_index(drop=True)
    result = target._deduplicate_stage1(df)
    pd.testing.assert_frame_equal(result, expected)


def test_sort_for_stage2_is_stable() -> None:
    df = pd.DataFrame(
        {
            "id": ["ID2", "ID1", "ID1"],
            "uniprot_id_primary": ["U1", "U1", "U1"],
            "target_chembl_id": ["C1", "C2", "C3"],
            "name": ["beta", "alpha", "gamma"],
        }
    )
    sorted_df = target._sort_for_stage2(df)
    assert list(sorted_df["id"]) == ["ID1", "ID1", "ID2"]
    assert list(sorted_df["name"])[:2] == ["alpha", "gamma"]


def test_dedup_stage2_distinct_three_columns() -> None:
    df = pd.DataFrame(
        {
            "id": ["ID1", "ID1", "ID2"],
            "uniprot_id_primary": ["U1", "U2", "U3"],
            "target_chembl_id": ["C1", "C1", "C2"],
            "name": ["alpha", "alpha", "beta"],
        }
    )
    result = target._deduplicate_stage2(df)
    assert len(result) == 2
    assert list(result["name"]) == ["alpha", "beta"]


def test_final_dedup_distinct_two_columns() -> None:
    df = pd.DataFrame(
        {
            "id": ["ID1", "ID1", "ID1"],
            "uniprot_id_primary": ["U1", "U2", "U3"],
            "target_chembl_id": ["C1", "C2", "C3"],
            "name": ["alpha", "alpha", "alpha"],
        }
    )
    result = target._deduplicate_final(df)
    assert len(result) == 1


def test_final_column_order(minimal_isoform_frame: pd.DataFrame) -> None:
    processed, _ = target._postprocess_frame(minimal_isoform_frame)
    assert list(processed.columns) == list(target.FINAL_COLUMN_ORDER)


def test_process_targets_output_name(tmp_path: Path, minimal_isoform_frame: pd.DataFrame) -> None:
    input_path = tmp_path / "output.target_20250101.csv"
    minimal_isoform_frame.to_csv(input_path, index=False)
    output_path = target.process_targets(str(input_path), verbose=False)
    assert output_path.name == "isoform.output.target_20250101.csv"


def test_process_targets_produces_expected_tokens(tmp_path: Path) -> None:
    input_path = tmp_path / "output.target_20250102.csv"
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "isoform_synonyms",
                "isoform_names",
                "isoform_ids",
                "uniprot_id_primary",
                "target_chembl_id",
            ]
        )
        writer.writerow(["PDE3A:alpha", "Alpha|A", "IsoA|IsoB", "Q99999", "CHEMBL999"])
    output_path = target.process_targets(str(input_path), verbose=False)
    df = pd.read_csv(output_path, dtype=str)
    assert set(df["name"]) == {"pde3a", "alpha", "3a", "a"}


def test_process_targets_deterministic(tmp_path: Path, snapshot_resource: Path) -> None:
    source = snapshot_resource / "targets_isoform_duplicates.csv"
    input_path = tmp_path / "output.target_20250103.csv"
    input_path.write_text(source.read_text(), encoding="utf-8")
    first = target.process_targets(str(input_path), verbose=False)
    second = target.process_targets(str(input_path), verbose=False)
    assert first.read_text(encoding="utf-8") == second.read_text(encoding="utf-8")


def test_make_triples_mismatched_lengths(snapshot_resource: Path) -> None:
    path = snapshot_resource / "targets_isoform_mismatched.csv"
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    triples = frame.apply(
        lambda row: target._make_triples(
            target._split_pipes(row["isoform_names"]),
            target._split_pipes(row["isoform_ids"]),
            target._split_pipes(row["isoform_synonyms"]),
        ),
        axis=1,
    )[0]
    assert len(triples) == 3
    assert triples[-1] == {"name": "Gamma", "id": None, "synonym": None}
