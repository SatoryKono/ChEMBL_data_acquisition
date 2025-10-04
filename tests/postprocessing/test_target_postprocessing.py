from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.config import Config
from library.common.csv_utils import write_csv_deterministic
from library.pipelines.target import postprocessing as target_postprocessing
from library.postprocessing import target
from library.schemas.targets import CELLULARITY_COLUMN_NAME, TARGETS_COLUMN_ORDER


@pytest.mark.unit
def test_split_pipes__removes_empty_segments() -> None:
    assert target._split_pipes("") == []
    assert target._split_pipes("a|b") == ["a", "b"]
    assert target._split_pipes(" a | | b ") == ["a", "b"]


@pytest.mark.unit
def test_matches_expected_input_name__accepts_suffixes() -> None:
    assert target._matches_expected_input_name("targets_20251002_chembl.csv")
    assert target._matches_expected_input_name("targets_20251002_chembl_normalized.csv")


@pytest.mark.unit
def test_make_triples__pads_shorter_sequences() -> None:
    triples = target._make_triples(
        names=["n1"],
        ids=["id1", "id2"],
        synonyms=["s1", "s2", "s3"],
    )
    assert triples == [
        {"name": "n1", "id": "id1", "synonym": "s1"},
        {"name": None, "id": "id2", "synonym": "s2"},
        {"name": None, "id": None, "synonym": "s3"},
    ]


@pytest.mark.unit
def test_lowercase_normalization__only_applies_to_names_and_synonyms(
    snapshot_resource: Path,
) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    ids = transform.result["id"].tolist()
    assert "ID_UP" in ids
    assert "id_low" in ids
    assert all(name == name.lower() for name in transform.result["name"])


@pytest.mark.unit
def test_tokenize_synonym__pde3a_alpha_variants() -> None:
    assert target._tokenize_synonym("PDE3A:alpha") == ["pde3a", "3a", "alpha"]


@pytest.mark.unit
def test_filter_names__removes_placeholder_entries(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    values = set(transform.result["name"].tolist())
    assert "" not in values
    assert "n/a" not in values
    assert "none" not in values


@pytest.mark.unit
def test_synonym_expansion__does_not_modify_ids(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    theta_row = transform.result.loc[transform.result["name"] == "theta"].iloc[0]
    assert theta_row["id"] == "ID_UP"


@pytest.mark.unit
def test_union_stage1__matches_distinct_on_four_columns(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    manual = transform.combined.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"]
    ).reset_index(drop=True)
    pdt.assert_frame_equal(transform.dedup_stage1.reset_index(drop=True), manual)


@pytest.mark.unit
def test_sorting_is_stable__uses_mergesort(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    manual = transform.dedup_stage1.sort_values(
        by=["uniprot_id_primary", "id"],
        kind="mergesort",
        na_position="first",
    ).reset_index(drop=True)
    pdt.assert_frame_equal(transform.sorted_stage, manual)


@pytest.mark.unit
def test_dedup_stage2__matches_distinct_on_three_columns(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    manual = transform.sorted_stage.drop_duplicates(
        subset=["id", "target_chembl_id", "name"], keep="first"
    ).reset_index(drop=True)
    pdt.assert_frame_equal(transform.dedup_stage2.reset_index(drop=True), manual)


@pytest.mark.unit
def test_final_distinct__drops_duplicates_on_two_columns(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    manual = transform.dedup_stage2.drop_duplicates(
        subset=["id", "name"], keep="first"
    ).reset_index(drop=True)
    pdt.assert_frame_equal(transform.result, manual)


@pytest.mark.unit
def test_result_column_order__matches_spec(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )
    transform = target._transform(frame)
    assert list(transform.result.columns) == [
        "id",
        "uniprot_id_primary",
        "target_chembl_id",
        "name",
    ]


@pytest.mark.unit
def test_process_targets__writes_expected_snapshot(tmp_path: Path, snapshot_resource: Path) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    expected = snapshot_resource / "target_isoform_minimal_expected.csv"
    working = tmp_path / "output.target_20250101.csv"
    shutil.copy(source, working)

    output_path = target.process_targets(str(working), verbose=False)

    result = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result, expected_frame)


@pytest.mark.unit
def test_output_filename__uses_isoform_prefix(tmp_path: Path, snapshot_resource: Path) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    working = tmp_path / "output.target_20251231.csv"
    shutil.copy(source, working)

    output_path = target.process_targets(str(working), verbose=False)
    assert output_path.name == "isoform.output.target_20251231.csv"


@pytest.mark.unit
def test_process_targets__auto_discovers_latest_input(tmp_path: Path, snapshot_resource: Path) -> None:
    first = tmp_path / "older.csv"
    first.write_text("dummy", encoding="utf-8")
    target_path = tmp_path / "output.target_20240101.csv"
    shutil.copy(snapshot_resource / "target_isoform_minimal.csv", target_path)
    newer = tmp_path / "output.target_20250101.csv"
    shutil.copy(snapshot_resource / "target_isoform_minimal.csv", newer)

    # Rewire search directory to the temporary folder.
    original_search_dir = target._DEFAULT_SEARCH_DIR
    target._DEFAULT_SEARCH_DIR = tmp_path
    try:
        output_path = target.process_targets(verbose=False)
    finally:
        target._DEFAULT_SEARCH_DIR = original_search_dir

    assert output_path.name == "isoform.output.target_20250101.csv"


@pytest.mark.unit
def test_process_targets__accepts_normalized_targets_file(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    expected = snapshot_resource / "target_isoform_minimal_expected.csv"
    working = tmp_path / "targets_20250101_normalized.csv"
    shutil.copy(source, working)

    output_path = target.process_targets(str(working), verbose=False)

    assert output_path.name == "isoform.targets_20250101_normalized.csv"
    result = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result, expected_frame)


@pytest.mark.unit
def test_finalise_targets__exports_power_query_column(
    tmp_path: Path, cfg: Config, snapshot_resource: Path
) -> None:
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_T1",
                "uniprotkb_Id": "P12345",
                "genus": "Homo",
                "lineage_superkingdom": "Eukaryota",
                "lineage_phylum": "Chordata",
                "lineage_class": "Mammalia",
                "superkingdom": "Eukaryota",
                "pref_name": "Target One",
                "type": "Legacy",
            },
            {
                "target_chembl_id": "CHEMBL_T2",
                "uniprotkb_Id": "Q99999",
                "genus": "Candida",
                "lineage_superkingdom": "Eukaryota",
                "lineage_phylum": "Ascomycota",
                "lineage_class": "Saccharomycetes",
                "superkingdom": "Eukaryota",
                "pref_name": "Target Two",
                "type": "Legacy",
            },
        ]
    )

    result = target_postprocessing.finalise_targets(frame)
    assert CELLULARITY_COLUMN_NAME in result.columns

    output_path = tmp_path / "targets_final.csv"
    write_csv_deterministic(
        result,
        output_path,
        col_order=TARGETS_COLUMN_ORDER,
        key_cols=["target_chembl_id"],
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
    )

    expected_path = snapshot_resource / "target_final_expected.csv"
    actual_bytes = output_path.read_bytes()
    expected_bytes = expected_path.read_bytes()
    assert actual_bytes == expected_bytes

    actual_frame = pd.read_csv(output_path, dtype=str, encoding=cfg.io.csv_encoding)
    expected_frame = pd.read_csv(
        expected_path, dtype=str, encoding=cfg.io.csv_encoding
    )
    pdt.assert_frame_equal(actual_frame, expected_frame)


@pytest.mark.unit
def test_process_targets__accepts_targets_file_with_suffix(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    expected = snapshot_resource / "target_isoform_minimal_expected.csv"
    working = tmp_path / "targets_20250101_chembl.csv"
    shutil.copy(source, working)

    output_path = target.process_targets(str(working), verbose=False)

    assert output_path.name == "isoform.targets_20250101_chembl.csv"
    result = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result, expected_frame)


@pytest.mark.unit
def test_process_targets__auto_discovers_latest_input_targets_prefix(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    first = tmp_path / "targets_20240101_normalized.csv"
    shutil.copy(snapshot_resource / "target_isoform_minimal.csv", first)
    newer = tmp_path / "targets_20250101_normalized.csv"
    shutil.copy(snapshot_resource / "target_isoform_minimal.csv", newer)

    original_search_dir = target._DEFAULT_SEARCH_DIR
    target._DEFAULT_SEARCH_DIR = tmp_path
    try:
        output_path = target.process_targets(verbose=False)
    finally:
        target._DEFAULT_SEARCH_DIR = original_search_dir

    assert output_path.name == "isoform.targets_20250101_normalized.csv"
