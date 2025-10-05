from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import (
    explain_exact_data_citation,
    explain_high_citation_rate,
    explain_higly_correlated_assay,
    explain_multmol_assay,
    explain_original_activity_flags,
    explain_review,
    explain_rounded_data_citation,
    explain_shuffled_assay,
    explain_unknown_chirality,
    process_activity_extended,
)
from library.postprocessing.activity_extended import (
    _ensure_required_input_columns,
    _rename_columns,
    _transform_activity_frame,
)

pytestmark = pytest.mark.postprocessing


@pytest.fixture()
def activity_flags_resources(snapshot_resource: Path) -> Path:
    return snapshot_resource / "activity_extended"


@pytest.fixture()
def activity_flags_exports(activity_flags_resources: Path, tmp_path: Path) -> Path:
    exports_src = activity_flags_resources / "exports"
    exports_dst = tmp_path / "exports"
    exports_dst.mkdir()
    for item in exports_src.iterdir():
        exports_dst.joinpath(item.name).write_bytes(item.read_bytes())
    dict_src = activity_flags_resources / "dictionary"
    dict_dst = tmp_path / "dictionary"
    dict_dst.mkdir()
    for subdir in dict_src.iterdir():
        target_dir = dict_dst / subdir.name
        target_dir.mkdir()
        for file in subdir.iterdir():
            target_dir.joinpath(file.name).write_bytes(file.read_bytes())
    return tmp_path


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_explain_unknown_chirality__cases() -> None:
    df = pd.DataFrame({"nstereo": [1, 2, pd.NA, 0]})
    explained = explain_unknown_chirality(df)
    assert list(explained["unknown_chirality"]) == [False, True, True, True]
    for message in explained["__explain_unknown_chirality"]:
        assert "unknown_chirality" in message


def test_explain_multmol_assay__duplicate_trigger() -> None:
    df = pd.DataFrame(
        {
            "salt_chembl_id": ["S1", "S1", "S2"],
            "molecule_chembl_id": ["M1", "M1", "M2"],
            "target_chembl_id": ["T1", "T1", "T2"],
            "assay_chembl_id": ["A1", "A1", "A2"],
            "standard_type": ["IC50", "IC50", "Ki"],
            "unknown_chirality": [False, False, True],
            "multmol_assay": [pd.NA, 0, 0],
        }
    )
    explained = explain_multmol_assay(df)
    assert list(explained["multmol_assay"]) == [True, True, False]
    assert "Count=2" in explained.loc[0, "__explain_multmol_assay"]
    assert "duplicate rule not triggered" in explained.loc[2, "__explain_multmol_assay"]


def test_boolean_flag_explainers__source_priority() -> None:
    df = pd.DataFrame(
        {
            "exact_cited_activity": [1, pd.NA],
            "exact_data_citation": [0, 1],
            "higly_correlated_cit": [0, pd.NA],
            "higly_correlated_assay": [1, 0],
            "shuffled_cit": [pd.NA, 1],
            "shuffled_assay": [0, 0],
            "review_doc": [1, pd.NA],
            "review": [pd.NA, 0],
            "approx_cited_activity": [pd.NA, 1],
            "rounded_data_citation": [0, pd.NA],
        }
    )
    exact_df = explain_exact_data_citation(df)
    assert list(exact_df["exact_data_citation"]) == [True, True]
    higly_df = explain_higly_correlated_assay(df)
    assert list(higly_df["higly_correlated_assay"]) == [False, False]
    assert list(higly_df["highly_correlated_assay"]) == [False, False]
    shuffled_df = explain_shuffled_assay(df)
    assert list(shuffled_df["shuffled_assay"]) == [False, True]
    review_df = explain_review(df)
    assert list(review_df["review"]) == [True, False]
    rounded_df = explain_rounded_data_citation(df)
    assert list(rounded_df["rounded_data_citation"]) == [False, True]


def test_explain_high_citation_rate__thresholds(activity_flags_resources: Path) -> None:
    raw = pd.read_csv(activity_flags_resources / "exports" / "output.activity_20251005.csv")
    renamed = _rename_columns(raw)
    explained = explain_high_citation_rate(
        renamed,
        dictionary_root=activity_flags_resources / "dictionary",
    )
    by_doc = explained.join(renamed["document_chembl_id"]).groupby("document_chembl_id")
    doc_flags = {
        doc: group["high_citation_rate"].unique().tolist()
        for doc, group in by_doc
    }
    assert doc_flags["DOC-1"] == [True]
    assert doc_flags["DOC-2"] == [False]
    assert doc_flags["DOC-3"] == [False]
    assert doc_flags["DOC-4"] == [False]


def test_explain_original_activity_flags__strings() -> None:
    df = pd.DataFrame(
        {
            "original_activity_approx": ["approx", pd.NA],
            "original_activity_exact": [pd.NA, "exact"],
        }
    )
    explained = explain_original_activity_flags(df)
    assert list(explained["original_activity_approx"]) == ["approx", pd.NA]
    assert list(explained["original_activity_exact"]) == [pd.NA, "exact"]
    for message in explained["__explain_original_activity_exact"]:
        assert "original_activity_exact" in message


def test_activity_flags_output_schema__dtypes(activity_flags_resources: Path) -> None:
    frame = pd.read_csv(activity_flags_resources / "exports" / "output.activity_20251005.csv")
    ensured, _ = _ensure_required_input_columns(frame)
    transformed = _transform_activity_frame(
        ensured,
        dictionary_root=activity_flags_resources / "dictionary",
        targets_override=None,
    )
    assert transformed["highly_correlated_assay"].dtype == transformed["higly_correlated_assay"].dtype
    for column in (
        "unknown_chirality",
        "multmol_assay",
        "exact_data_citation",
        "higly_correlated_assay",
        "highly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
        "high_citation_rate",
    ):
        assert str(transformed[column].dtype) == "boolean"
    assert str(transformed["original_activity_approx"].dtype) == "object"
    assert str(transformed["original_activity_exact"].dtype) == "object"


EXPECTED_ACTIVITY_HASH = "7ed915aec9d6cbf75585722750474d5b7444f4070ad1721f972d3490e237671b"


def test_process_activity_extended__golden_snapshot(activity_flags_exports: Path) -> None:
    exports = activity_flags_exports / "exports"
    dictionary = activity_flags_exports / "dictionary"
    input_path = exports / "output.activity_20251005.csv"
    output_path = process_activity_extended(
        input_path=input_path,
        dictionary_dir=dictionary,
    )
    expected = pd.read_csv(Path("tests/golden/extended.output.activity_20251005.csv"))
    result = pd.read_csv(output_path)
    pd.testing.assert_frame_equal(result[expected.columns], expected)
    digest = _file_sha256(output_path)
    if EXPECTED_ACTIVITY_HASH:
        assert digest == EXPECTED_ACTIVITY_HASH
