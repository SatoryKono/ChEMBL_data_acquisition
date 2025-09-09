from __future__ import annotations

import pandas as pd
import sys
from pathlib import Path

import pytest

from library.input_initialisation_library import (
    _ensure_openpyxl,
    TableDict,
    append_entities,
    build_combined_tables,
    unify_dtypes,
    save_tables,
    process_activity_table,
)
import library.input_initialisation_library as lib


def test_unify_dtypes_basic() -> None:
    df = pd.DataFrame(
        {
            "assay_chembl_id": ["A1"],
            "isoform": ["1"],
            "shuffled_cit": ["1"],
            "mw_freebase": ["10.5"],
            "publication_date": ["2020-01-01"],
        }
    )
    res = unify_dtypes(df)
    assert res["assay_chembl_id"].dtype == "string"
    assert str(res["isoform"].dtype) == "Int64"
    assert res["shuffled_cit"].dtype == "boolean"
    assert res["mw_freebase"].dtype == "float64"
    assert str(res["publication_date"].dtype).startswith("datetime64[")


def test_append_entities_deduplication() -> None:
    df_a = pd.DataFrame({"id": [1, 2], "val": ["a", "b"]})
    df_b = pd.DataFrame({"id": [2, 3], "val": ["b", "c"]})
    res = append_entities(df_a, df_b)
    assert len(res) == 3
    assert sorted(res["id"].tolist()) == [1, 2, 3]


def test_build_combined_tables_drops_activity_cols() -> None:
    same: TableDict = {
        "activity": pd.DataFrame({"id": [1], "Column1": ["a"]}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame({"id": [2], "Column1": ["b"]}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs": pd.DataFrame(),
    }
    combined = build_combined_tables(same, all_)
    assert "Column1" not in combined["activity"].columns
    assert len(combined["activity"]) == 2
    assert "pairs_same_document" in combined
    assert "pairs" in combined


def test_build_combined_tables_pairs_no_merge() -> None:
    """Pairs from different sources should not be combined."""
    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame({"id": range(3)}),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame(),
        "pairs": pd.DataFrame({"id": range(5)}),
    }
    combined = build_combined_tables(same, all_)
    assert len(combined["pairs_same_document"]) == 3
    assert len(combined["pairs"]) == 5


def test_build_combined_tables_handles_status(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``build_combined_tables`` should populate activity before status logic."""
    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [1]}),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [2]}),
        "pairs": pd.DataFrame(),
    }

    # Minimal status.csv to satisfy loader
    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n" "ok,null,null,0,0\n"
    )

    # Simplify heavy processing steps
    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir: df)
    monkeypatch.setattr(
        lib,
        "initialize_activity_status",
        lambda df, _api: df.assign(**{"Filtered.init": "ok"}),
    )

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    assert len(combined["activity"]) == 2
    assert "Filtered.init" in combined["activity"].columns


def test_build_combined_tables_initializes_pair_tables(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Pair tables should receive Filtered columns via ``initialize_pairs``."""
    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [1]}),
        "pairs_same_document": pd.DataFrame({"activity_id1": [1], "activity_id2": [2]}),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [2]}),
        "pairs": pd.DataFrame({"activity_id1": [2], "activity_id2": [1]}),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n"
        "good,null,null,0,0\n"
        "bad,null,null,1,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir: df)

    def fake_init(df: pd.DataFrame, _api: object) -> pd.DataFrame:
        mapping = {1: "good", 2: "bad"}
        return df.assign(**{"Filtered.init": df["activity_id"].map(mapping)})

    monkeypatch.setattr(lib, "initialize_activity_status", fake_init)
    monkeypatch.setattr(lib, "aggregate_activity", lambda *_: {})

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    pairs = combined["pairs"].iloc[0]
    assert [pairs["Filtered1"], pairs["Filtered2"], pairs["Filtered"]] == [
        "bad",
        "good",
        "good",
    ]
    pairs_same = combined["pairs_same_document"].iloc[0]
    assert [
        pairs_same["Filtered1"],
        pairs_same["Filtered2"],
        pairs_same["Filtered"],
    ] == [
        "good",
        "bad",
        "good",
    ]


def test_build_combined_tables_normalizes_pair_column_names(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Pair tables with varied activity ID column names should still initialise."""

    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [1]}),
        "pairs_same_document": pd.DataFrame({"Activity_ID1": [1], "Activity_ID2": [2]}),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_id": [2]}),
        "pairs": pd.DataFrame({"ACTIVITY_ID_1": [2], "ACTIVITY_ID_2": [1]}),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n" "good,null,null,0,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir: df)

    def fake_init(df: pd.DataFrame, _api: object) -> pd.DataFrame:
        mapping = {1: "good", 2: "good"}
        return df.assign(**{"Filtered.init": df["activity_id"].map(mapping)})

    monkeypatch.setattr(lib, "initialize_activity_status", fake_init)
    monkeypatch.setattr(lib, "aggregate_activity", lambda *_: {})

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    for key in ("pairs", "pairs_same_document"):
        assert {
            "activity_chembl_id1",
            "activity_chembl_id2",
            "Filtered1",
            "Filtered2",
            "Filtered",
        }.issubset(combined[key].columns)


def test_save_tables_writes_files(tmp_path: Path) -> None:
    tables: TableDict = {
        "activity": pd.DataFrame({"id": [1]}),
        "assay": pd.DataFrame({"id": [2]}),
        "document": pd.DataFrame({"id": [3]}),
        "target": pd.DataFrame({"id": [4]}),
        "testitem": pd.DataFrame({"id": [5]}),
        "pairs_same_document": pd.DataFrame({"id": [1]}),
        "pairs": pd.DataFrame({"id": [1]}),
    }
    paths = save_tables(tables, tmp_path)
    for entity, path in paths.items():
        assert path.exists(), f"missing {entity} file"
        df = pd.read_csv(path)
        assert not df.empty


def test_ensure_openpyxl_version(monkeypatch: pytest.MonkeyPatch) -> None:
    module = type("M", (), {"__version__": "3.0.0"})()
    monkeypatch.setitem(sys.modules, "openpyxl", module)
    with pytest.raises(RuntimeError):
        _ensure_openpyxl()
    monkeypatch.delitem(sys.modules, "openpyxl")


def test_process_activity_table_basic(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "activity_chembl_id": [1, 2],
            "salt_chembl_id": ["S1", "S1"],
            "molecule_chembl_id": ["M1", "M1"],
            "target_chembl_id": ["T1", "T1"],
            "assay_chembl_id": ["A1", "A1"],
            "document_chembl_id": ["D1", "D1"],
            "bao_endpoint": ["e1", "e1"],
            "standard_type": ["t", "t"],
            "standard_value": [1.0, 2.0],
            "log_value": [0.1, 0.2],
            "bao_format": ["f1", "f1"],
            "compound_key": ["ck1", "ck1"],
            "compound_name": ["cn1", "cn1"],
            "multmol_assay": [pd.NA, pd.NA],
            "nstereo": [2, 1],
            "approx_cited_activity": [pd.NA, False],
            "exact_cited_activity": [True, False],
            "higly_correlated_cit": [pd.NA, False],
            "shuffled_cit": [pd.NA, False],
            "review_doc": [False, False],
            "original_activity_approx": [0.5, 0.6],
            "original_activity_exact": [0.5, 0.6],
        }
    )

    (tmp_path / "citation_fraction.csv").write_text(
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n2,1,x,0.05\n"
    )
    (tmp_path / "targets_type.csv").write_text(
        "chembl_id,type,IUPHAR_class,IUPHAR_subclass\nT1,Unicellular organism,ClassA,Multifunctional\n"
    )

    res = process_activity_table(df, tmp_path)
    expected_cols = [
        "activity_id",
        "saltform_id",
        "molecule_id",
        "target_id",
        "assay_id",
        "document_id",
        "bao_endpoint",
        "mesurement_type",
        "standard_value",
        "pA_value",
        "bao_format",
        "compound_key",
        "compound_name",
        "unknown_chirality",
        "multmol_assay",
        "exact_data_citation",
        "higly_correlated_assay",
        "shuffled_assay",
        "review",
        "rounded_data_citation",
        "high_citation_rate",
        "original_activity_approx",
        "original_activity_exact",
        "is_citation",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "unicellular_organism",
        "multifunctional_enzyme",
    ]
    assert list(res.columns) == expected_cols
    assert bool(res.loc[0, "is_citation"])
    assert not bool(res.loc[1, "is_citation"])
    assert res["high_citation_rate"].all()
    assert res["unicellular_organism"].all()

    assert res["multifunctional_enzyme"].all()


def test_process_activity_table_without_nstereo(tmp_path: Path) -> None:
    """Process activity table lacking ``nstereo`` column."""
    df = pd.DataFrame(
        {
            "activity_chembl_id": [1],
            "salt_chembl_id": ["S1"],
            "molecule_chembl_id": ["M1"],
            "target_chembl_id": ["T1"],
            "assay_chembl_id": ["A1"],
            "document_chembl_id": ["D1"],
            "bao_endpoint": ["e1"],
            "standard_type": ["t"],
            "standard_value": [1.0],
            "log_value": [0.1],
            "bao_format": ["f1"],
            "compound_key": ["ck1"],
            "compound_name": ["cn1"],
            "multmol_assay": [pd.NA],
            "approx_cited_activity": [pd.NA],
            "exact_cited_activity": [True],
            "higly_correlated_cit": [pd.NA],
            "shuffled_cit": [pd.NA],
            "review_doc": [False],
            "original_activity_approx": [0.5],
            "original_activity_exact": [0.5],
        }
    )

    (tmp_path / "citation_fraction.csv").write_text(
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n1,1,x,0.05\n"
    )
    (tmp_path / "targets_type.csv").write_text(
        "chembl_id,type,IUPHAR_class,IUPHAR_subclass\nT1,Multicellular organism,,\n"
    )

    res = process_activity_table(df, tmp_path)
    assert "unknown_chirality" in res.columns
    assert res.loc[0, "unknown_chirality"]

    assert pd.isna(res.loc[0, "multifunctional_enzyme"])
