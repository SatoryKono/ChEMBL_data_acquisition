from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import target


def _copy_resource(src: Path, dst: Path) -> Path:
    dst.write_bytes(src.read_bytes())
    return dst


@pytest.mark.unit
def test_normalise_series__strips_sentinels_and_whitespace() -> None:
    series = pd.Series(["  Value  ", "n/a", None, "N.A", "n / a", "None"])

    result = target._normalise_series(series, lowercase=False)

    assert result.tolist() == ["Value", "", "", "", "", ""]


@pytest.mark.unit
def test_normalise_series__lowercases_when_requested() -> None:
    series = pd.Series(["Alpha", " Beta "])

    result = target._normalise_series(series, lowercase=True)

    assert result.tolist() == ["alpha", "beta"]


@pytest.mark.unit
def test_split_pipe__filters_empty_segments() -> None:
    assert target._split_pipe("foo | | bar|| baz ") == ["foo", "bar", "baz"]


@pytest.mark.unit
def test_split_synonyms__deduplicates_and_ignores_na_tokens() -> None:
    value = "Alpha:alpha:beta:NA:none:n/a"

    result = target._split_synonyms(value)

    assert result == ["Alpha", "alpha", "beta"]


@pytest.mark.unit
def test_syn_expand__splits_on_non_alphanumeric_boundaries() -> None:
    term = "Alpha-beta's isoform"

    result = target._syn_expand(term)

    assert result == ["alpha", "beta's", "isoform"]


@pytest.mark.unit
def test_syn_expand__returns_empty_for_blank_input() -> None:
    assert target._syn_expand("") == []


@pytest.mark.unit
def test_process_targets__writes_default_filename(tmp_path: Path, snapshot_resource: Path) -> None:
    source = _copy_resource(
        snapshot_resource / "target_postprocessing_input.csv",
        tmp_path / "targets.csv",
    )

    output_path = target.process_targets(source)

    expected_path = tmp_path / "isoform.output.targets.csv"
    assert output_path == expected_path
    assert output_path.exists()


@pytest.mark.unit
def test_process_targets__respects_output_directory(tmp_path: Path, snapshot_resource: Path) -> None:
    source = _copy_resource(
        snapshot_resource / "target_postprocessing_input.csv",
        tmp_path / "targets.csv",
    )
    destination = tmp_path / "out"

    output_path = target.process_targets(source, output_dir=destination)

    expected_path = destination / "isoform.output.targets.csv"
    assert output_path == expected_path
    assert output_path.exists()


@pytest.mark.unit
def test_process_targets__expands_isoform_entries_into_token_rows(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = _copy_resource(
        snapshot_resource / "target_postprocessing_input.csv",
        tmp_path / "targets.csv",
    )

    output_path = target.process_targets(source)
    result = pd.read_csv(output_path, keep_default_na=False)

    assert list(result.columns) == [
        "target_chembl_id",
        "isoform_id",
        "isoform_name",
        "term",
        "token",
    ]

    alpha_rows = result[(result["target_chembl_id"] == "CHEMBLT1") & (result["isoform_id"] == "ENSP000001")]
    assert set(alpha_rows[alpha_rows["term"] == "Alpha Isoform"]["token"]) == {"alpha", "isoform"}
    assert (alpha_rows[alpha_rows["term"] == "alpha chain"]["token"].tolist()) == ["alpha", "chain"]

    synonym_only = result[
        (result["target_chembl_id"] == "CHEMBLT3")
        & (result["isoform_id"] == "")
        & (result["term"] == "functional isoform")
    ]
    assert synonym_only["token"].tolist() == ["functional", "isoform"]


@pytest.mark.unit
def test_process_targets__skips_rows_without_names_or_synonyms(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = _copy_resource(
        snapshot_resource / "target_postprocessing_input.csv",
        tmp_path / "targets.csv",
    )

    output_path = target.process_targets(source)
    result = pd.read_csv(output_path, keep_default_na=False)

    assert "CHEMBLT4" not in result["target_chembl_id"].unique()


@pytest.mark.unit
def test_process_targets__deduplicates_duplicate_records(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBLD1",
                "isoform_ids": "ENSP000010",
                "isoform_names": "iso",
                "isoform_synonyms": "iso",
            },
            {
                "target_chembl_id": "CHEMBLD1",
                "isoform_ids": "ENSP000010",
                "isoform_names": "iso",
                "isoform_synonyms": "iso",
            },
        ]
    )
    source = tmp_path / "duplicates.csv"
    df.to_csv(source, index=False)

    output_path = target.process_targets(source)
    result = pd.read_csv(output_path, keep_default_na=False)

    assert len(result) == 1
    row = result.loc[0]
    assert row.to_dict() == {
        "target_chembl_id": "CHEMBLD1",
        "isoform_id": "ENSP000010",
        "isoform_name": "iso",
        "term": "iso",
        "token": "iso",
    }


@pytest.mark.unit
def test_process_targets__sorts_rows_by_identifiers(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBLA",
                "isoform_ids": "ENSP2|ENSP1",
                "isoform_names": "b|a",
                "isoform_synonyms": "b:b|a:a",
            }
        ]
    )
    source = tmp_path / "unsorted.csv"
    df.to_csv(source, index=False)

    output_path = target.process_targets(source)
    result = pd.read_csv(output_path, keep_default_na=False)
    sorted_result = result.sort_values(
        by=["target_chembl_id", "isoform_id", "term", "token"],
        kind="mergesort",
        ignore_index=True,
    )

    pd.testing.assert_frame_equal(result, sorted_result)


@pytest.mark.unit
def test_process_targets__reads_using_fallback_encoding(tmp_path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBLCYR",
                "isoform_ids": "ENSPRU",
                "isoform_names": "йод",
                "isoform_synonyms": "йод",
            }
        ]
    )
    source = tmp_path / "cyrillic.csv"
    df.to_csv(source, index=False, encoding="cp1251")

    output_path = target.process_targets(source)
    result = pd.read_csv(output_path, keep_default_na=False)

    assert result.empty
    assert list(result.columns) == [
        "target_chembl_id",
        "isoform_id",
        "isoform_name",
        "term",
        "token",
    ]

