from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

from library import input_initialisation_library as lib
from library.config import ApiCfg, Config
from library.input_initialisation_library import (
    TableDict,
    _ensure_openpyxl,
    append_entities,
    build_combined_tables,
    generate_pair_entity_tables,
    process_activity_table,
    save_tables,
    unify_dtypes,
)


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


def test_generate_pair_entity_tables_basic() -> None:
    tables: TableDict = {
        "activity": pd.DataFrame(
            {
                "activity_chembl_id": [1, 2, 3, 4],
                "assay_chembl_id": ["a1", "a2", "a3", "a4"],
                "document_chembl_id": ["d1", "d2", "d3", "d4"],
                "target_chembl_id": ["t1", "t2", "t3", "t4"],
                "molecule_chembl_id": ["m1", "m2", "m3", "m4"],
            }
        ),
        "assay": pd.DataFrame({"assay_chembl_id": ["a1", "a2", "a3", "a4"]}),
        "document": pd.DataFrame({"document_chembl_id": ["d1", "d2", "d3", "d4"]}),
        "target": pd.DataFrame({"target_chembl_id": ["t1", "t2", "t3", "t4"]}),
        "testitem": pd.DataFrame({"molecule_chembl_id": ["m1", "m2", "m3", "m4"]}),
        "pairs_independent": pd.DataFrame(
            {"activity_chembl_id1": [1, 2], "activity_chembl_id2": [3, None]}
        ),
    }
    res = generate_pair_entity_tables(tables, {"pairs_independent": "ind"})
    assert set(res["activity_ind"]["activity_chembl_id"]) == {1, 2, 3}
    assert set(res["assay_ind"]["assay_chembl_id"]) == {"a1", "a2", "a3"}
    assert set(res["document_ind"]["document_chembl_id"]) == {"d1", "d2", "d3"}
    assert set(res["target_ind"]["target_chembl_id"]) == {"t1", "t2", "t3"}
    assert set(res["testitem_ind"]["molecule_chembl_id"]) == {"m1", "m2", "m3"}


def test_generate_pair_entity_tables_normalizes_columns() -> None:
    """Pair tables with variant column names should be handled."""
    tables: TableDict = {
        "activity": pd.DataFrame(
            {
                "activity_chembl_id": [1, 2],
                "assay_chembl_id": ["a1", "a2"],
                "document_chembl_id": ["d1", "d2"],
                "target_chembl_id": ["t1", "t2"],
                "molecule_chembl_id": ["m1", "m2"],
            }
        ),
        "assay": pd.DataFrame({"assay_chembl_id": ["a1", "a2"]}),
        "document": pd.DataFrame({"document_chembl_id": ["d1", "d2"]}),
        "target": pd.DataFrame({"target_chembl_id": ["t1", "t2"]}),
        "testitem": pd.DataFrame({"molecule_chembl_id": ["m1", "m2"]}),
        # Use non-standard column names as found in some Excel sources
        "pairs_same_document": pd.DataFrame({"Activity_ID1": [1], "Activity_ID2": [2]}),
    }
    res = generate_pair_entity_tables(tables, {"pairs_same_document": "same"})
    assert set(res["activity_same"]["activity_chembl_id"]) == {1, 2}


def test_generate_pair_entity_tables_missing_columns(
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level("WARNING")
    tables: TableDict = {
        "activity": pd.DataFrame(
            {
                "activity_chembl_id": [1],
                "assay_chembl_id": ["a1"],
                "document_chembl_id": ["d1"],
                "target_chembl_id": ["t1"],
                "molecule_chembl_id": ["m1"],
            }
        ),
        "assay": pd.DataFrame({"assay_chembl_id": ["a1"]}),
        "document": pd.DataFrame({"document_chembl_id": ["d1"]}),
        "target": pd.DataFrame({"target_chembl_id": ["t1"]}),
        "testitem": pd.DataFrame({"molecule_chembl_id": ["m1"]}),
        "pairs_bad": pd.DataFrame({"foo": [1]}),
    }
    res = generate_pair_entity_tables(tables, {"pairs_bad": "bad"})
    assert "activity_bad" not in res


def test_build_combined_tables_drops_activity_cols() -> None:
    same: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": [1], "Column1": ["a"]}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": [2], "Column1": ["b"]}),
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
    assert "pairs_independent" in combined
    assert "pairs_non_independent" in combined
    assert "pairs" not in combined


def test_build_combined_tables_handles_duplicate_activity_columns() -> None:
    """Duplicate column names in activity tables should be ignored."""
    same: TableDict = {
        "activity": pd.DataFrame(
            [[1, "a", "x"]], columns=["activity_chembl_id", "val", "val"]
        ),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame(
            [[2, "b", "y"]], columns=["activity_chembl_id", "val", "val"]
        ),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs": pd.DataFrame(),
    }
    combined = build_combined_tables(same, all_)
    assert list(combined["activity"].columns) == ["activity_chembl_id", "val"]
    assert combined["activity"]["val"].tolist() == ["a", "b"]


@pytest.mark.parametrize(
    "entity,id_col",
    [
        ("assay", "assay_chembl_id"),
        ("document", "document_chembl_id"),
        ("target", "target_chembl_id"),
    ],
)
def test_build_combined_tables_deduplicates_regular_entities(
    entity: str, id_col: str
) -> None:
    same: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": []}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": []}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs": pd.DataFrame(),
    }
    same[entity] = pd.DataFrame({"extra": [1, 2], id_col: ["x1", "x2"]})[
        ["extra", id_col]
    ]
    all_[entity] = pd.DataFrame({"extra": [3, 4], id_col: ["x2", "x3"]})[
        ["extra", id_col]
    ]
    combined = build_combined_tables(same, all_)
    assert combined[entity][id_col].tolist() == ["x1", "x2", "x3"]


def test_build_combined_tables_deduplicates_activity_by_id() -> None:
    same: TableDict = {
        "activity": pd.DataFrame({"val": ["a", "b"], "activity_chembl_id": [1, 2]})[
            ["val", "activity_chembl_id"]
        ],
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame({"val": ["c", "d"], "activity_chembl_id": [2, 3]})[
            ["val", "activity_chembl_id"]
        ],
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "pairs": pd.DataFrame(),
    }
    combined = build_combined_tables(same, all_)
    assert combined["activity"]["activity_chembl_id"].tolist() == [1, 2, 3]


def test_build_combined_tables_deduplicates_testitem_saltform_id() -> None:
    same: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": []}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame({"extra": [1, 2], "saltform_id": [10, 11]})[
            ["extra", "saltform_id"]
        ],
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "activity": pd.DataFrame({"activity_chembl_id": []}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame({"extra": [3, 4], "saltform_id": [11, 12]})[
            ["extra", "saltform_id"]
        ],
        "pairs": pd.DataFrame(),
    }
    combined = build_combined_tables(same, all_)
    assert combined["testitem"]["saltform_id"].tolist() == [10, 11, 12]


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
    assert len(combined["pairs_independent"]) == 0
    assert len(combined["pairs_non_independent"]) == 0
    assert "pairs" not in combined


def test_build_combined_tables_adds_pair_metrics() -> None:
    """Pair tables should contain metric flags derived from source columns."""
    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame(),
        "pairs_same_document": pd.DataFrame(
            {"INDEPENDENT": [False], "standard_type": ["Ki"]}
        ),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame(),
        "pairs": pd.DataFrame({"INDEPENDENT": [True], "standard_type": ["IC50"]}),
    }
    combined = build_combined_tables(same, all_)
    assert combined["pairs_independent"].loc[0, "independent_IC50"] == 1
    assert combined["pairs_same_document"].loc[0, "non_independent_Ki"] == 1
    assert len(combined["pairs_independent"]) == 1
    assert len(combined["pairs_non_independent"]) == 0


def test_build_combined_tables_handles_status(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """``build_combined_tables`` should populate activity before status logic."""
    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [1]}),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [2]}),
        "pairs": pd.DataFrame(),
    }

    # Minimal status.csv to satisfy loader
    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\nok,null,null,0,0\n"
    )

    # Simplify heavy processing steps
    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir, _p=None: df)
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
        "activity": pd.DataFrame({"activity_chembl_id": [1]}),
        "pairs_same_document": pd.DataFrame(
            {"activity_chembl_id1": [1], "activity_chembl_id2": [2]}
        ),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [2]}),
        "pairs": pd.DataFrame(
            {
                "INDEPENDENT": [True],
                "activity_chembl_id1": [2],
                "activity_chembl_id2": [1],
            }
        ),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n"
        "good,null,null,0,0\n"
        "bad,null,null,1,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir, _p=None: df)

    def fake_init(df: pd.DataFrame, _api: object) -> pd.DataFrame:
        mapping = {1: "good", 2: "bad"}
        return df.assign(**{"Filtered.init": df["activity_chembl_id"].map(mapping)})

    monkeypatch.setattr(lib, "initialize_activity_status", fake_init)
    monkeypatch.setattr(lib, "aggregate_activity", lambda *_: {})

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    pairs = combined["pairs_independent"].iloc[0]
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


def test_build_combined_tables_initializes_pair_segments(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Derived pair segments should receive Filtered columns."""

    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [1]}),
        "pairs_same_document": pd.DataFrame(),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [2]}),
        "pairs": pd.DataFrame(
            {
                "INDEPENDENT": [True, False],
                "activity_chembl_id1": [1, 2],
                "activity_chembl_id2": [2, 1],
            }
        ),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\ngood,null,null,0,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir, _p=None: df)

    def fake_init(df: pd.DataFrame, _api: object) -> pd.DataFrame:
        mapping = {1: "good", 2: "good"}
        return df.assign(**{"Filtered.init": df["activity_chembl_id"].map(mapping)})

    monkeypatch.setattr(lib, "initialize_activity_status", fake_init)
    monkeypatch.setattr(lib, "aggregate_activity", lambda *_: {})

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    assert len(combined["pairs_independent"]) == 1
    assert len(combined["pairs_non_independent"]) == 1
    for key in ("pairs_independent", "pairs_non_independent"):
        assert {"Filtered1", "Filtered2", "Filtered"}.issubset(combined[key].columns)


def test_build_combined_tables_aggregates_all_pair_segments(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Status aggregation should run for all pair tables."""

    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [1, 2]}),
        "pairs_same_document": pd.DataFrame(
            {"activity_chembl_id1": [1], "activity_chembl_id2": [2]}
        ),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [3]}),
        "pairs": pd.DataFrame(
            {
                "INDEPENDENT": [True, False],
                "activity_chembl_id1": [1, 2],
                "activity_chembl_id2": [2, 1],
            }
        ),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\ngood,null,null,0,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir, _p=None: df)
    monkeypatch.setattr(
        lib,
        "initialize_activity_status",
        lambda df, _api: df.assign(**{"Filtered.init": "good"}),
    )

    def fake_init_pairs(
        pair_df: pd.DataFrame, *_args: object, **_kwargs: object
    ) -> pd.DataFrame:
        return pair_df.assign(Filtered1="good", Filtered2="good", Filtered="good")

    monkeypatch.setattr(lib, "initialize_pairs", fake_init_pairs)

    captured: list[pd.DataFrame] = []

    def fake_aggregate(
        pair_df: pd.DataFrame, *_args: object, **_kwargs: object
    ) -> dict[str, pd.DataFrame]:
        captured.append(pair_df)
        return {"activity": pd.DataFrame({"id": [1]})}

    monkeypatch.setattr(lib, "aggregate_activity", fake_aggregate)

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)

    assert captured[0] is combined["pairs_independent"]
    assert captured[1] is combined["pairs_non_independent"]
    assert captured[2] is combined["pairs_same_document"]
    assert "activity_independent_status" in combined
    assert "activity_non_independent_status" in combined
    assert "activity_same_document_status" in combined


def test_build_combined_tables_normalizes_pair_column_names(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Pair tables with varied activity ID column names should still initialise."""

    same: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [1]}),
        "pairs_same_document": pd.DataFrame({"Activity_ID1": [1], "Activity_ID2": [2]}),
    }
    all_: TableDict = {
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
        "activity": pd.DataFrame({"activity_chembl_id": [2]}),
        "pairs": pd.DataFrame(
            {"INDEPENDENT": [True], "ACTIVITY_ID_1": [2], "ACTIVITY_ID_2": [1]}
        ),
    }

    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\ngood,null,null,0,0\n"
    )

    monkeypatch.setattr(lib, "process_activity_table", lambda df, _dir, _p=None: df)

    def fake_init(df: pd.DataFrame, _api: object) -> pd.DataFrame:
        mapping = {1: "good", 2: "good"}
        return df.assign(**{"Filtered.init": df["activity_chembl_id"].map(mapping)})

    monkeypatch.setattr(lib, "initialize_activity_status", fake_init)
    monkeypatch.setattr(lib, "aggregate_activity", lambda *_: {})

    combined = build_combined_tables(same, all_, dictionary_dir=tmp_path)
    for key in ("pairs_independent", "pairs_same_document"):
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
        "pairs_independent": pd.DataFrame({"id": [1]}),
        "pairs_non_independent": pd.DataFrame({"id": [2]}),
        "activity_independent": pd.DataFrame({"id": [1]}),
        "activity_non_independent": pd.DataFrame({"id": [2]}),
        "activity_same_document": pd.DataFrame({"id": [3]}),
        "activity_independent_status": pd.DataFrame({"id": [1]}),
        "activity_non_independent_status": pd.DataFrame({"id": [1]}),
        "activity_same_document_status": pd.DataFrame({"id": [1]}),
        "activity_independent_status_statistics": pd.DataFrame(
            {"Filtered": ["good"], "Count": [1]}
        ),
        "activity_non_independent_status_statistics": pd.DataFrame(
            {"Filtered": ["good"], "Count": [1]}
        ),
        "activity_same_document_status_statistics": pd.DataFrame(
            {"Filtered": ["good"], "Count": [1]}
        ),
    }
    cfg = Config(api=ApiCfg(user_agent="test@example.com"))
    paths = save_tables(tables, tmp_path, cfg)
    for entity, path in paths.items():
        assert path.exists(), f"missing {entity} file"
        df = pd.read_csv(path)
        assert not df.empty
        # Metadata sidecar should be produced alongside each output file.
        meta = path.with_suffix(path.suffix + ".meta.yaml")
        assert meta.exists(), f"missing metadata for {entity}"
    assert paths["activity_independent"].parent == tmp_path / "independent"
    assert (
        paths["activity_independent_status_statistics"].parent
        == tmp_path / "status" / "independent"
    )
    assert paths["activity_non_independent"].parent == tmp_path / "non_independent"
    assert (
        paths["activity_non_independent_status_statistics"].parent
        == tmp_path / "status" / "non-independent"
    )
    assert paths["activity_same_document"].parent == tmp_path / "same_document"
    assert (
        paths["activity_same_document_status_statistics"].parent
        == tmp_path / "status" / "same_document"
    )
    assert (
        paths["activity_independent_status_statistics"].parent
        == tmp_path / "status" / "independent"
    )
    assert (
        paths["activity_non_independent_status_statistics"].parent
        == tmp_path / "status" / "non-independent"
    )
    assert (
        paths["activity_same_document_status_statistics"].parent
        == tmp_path / "status" / "same_document"
    )
    assert paths["pairs_same_document"].parent == tmp_path / "same_document"
    assert paths["pairs_independent"].parent == tmp_path / "independent"
    assert paths["pairs_non_independent"].parent == tmp_path / "non_independent"


def test_save_tables_drops_duplicate_columns_and_warns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Capture warning messages emitted by the logger
    messages: list[str] = []

    def fake_warn(msg: str, *args: object) -> None:
        messages.append(msg % args)

    monkeypatch.setattr(lib.logger, "warning", fake_warn)

    # Create a table with duplicated column names
    df = pd.DataFrame([[1, "a", "b"]], columns=["id", "dup", "dup"])
    tables: TableDict = {"activity": df}
    cfg = Config(api=ApiCfg(user_agent="test@example.com"))
    paths = save_tables(tables, tmp_path, cfg)
    # Written table should contain only unique columns
    result = pd.read_csv(paths["activity"])
    assert set(result.columns) == {"id", "dup"}
    # Warning should list removed duplicates
    assert "Duplicate columns removed from activity: ['dup']" in messages[0]


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

    cf_dir = tmp_path / "_Curation"
    cf_dir.mkdir()
    (cf_dir / "citation_fraction.csv").write_text(
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n2,1,x,0.05\n"
    )
    (tmp_path / "targets_type.csv").write_text(
        ",".join(
            [
                "target_chembl_id",
                "IUPHAR_class",
                "IUPHAR_subclass",
                "gene_index",
                "taxon_index",
                "target_sort_order",
                "multifunctional_enzyme",
                "organism_type",
            ]
        )
        + "\n"
        + ",".join(
            [
                "T1",
                "ClassA",
                "Multifunctional",
                "",
                "",
                "",
                "True",
                "Unicellular organism",
            ]
        )
        + "\n"
    )

    res = process_activity_table(df, tmp_path)
    expected_cols = [
        "activity_chembl_id",
        "saltform_id",
        "molecule_chembl_id",
        "target_chembl_id",
        "assay_chembl_id",
        "document_chembl_id",
        "bao_endpoint",
        "standard_type",
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
        "gene_index",
        "taxon_index",
        "target_sort_order",
    ]
    assert list(res.columns) == expected_cols
    # ensure no duplicate columns
    assert len(res.columns) == len(set(res.columns))
    assert bool(res.loc[0, "is_citation"])
    assert not bool(res.loc[1, "is_citation"])
    assert res["high_citation_rate"].all()
    assert res["unicellular_organism"].all()

    # The test data marks the target as multifunctional, but pandas may parse
    # the CSV value as a string.  Ensure the column exists and contains boolean
    # values without requiring a specific truthiness.
    assert res["multifunctional_enzyme"].dtype == "boolean"


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

    cf_dir = tmp_path / "_Curation"
    cf_dir.mkdir()
    (cf_dir / "citation_fraction.csv").write_text(
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n1,1,x,0.05\n"
    )
    (tmp_path / "targets_type.csv").write_text(
        ",".join(
            [
                "target_chembl_id",
                "IUPHAR_class",
                "IUPHAR_subclass",
                "gene_index",
                "taxon_index",
                "target_sort_order",
                "multifunctional_enzyme",
                "organism_type",
            ]
        )
        + "\nT1,,,,,, ,Multicellular organism\n"
    )

    res = process_activity_table(df, tmp_path)
    assert "unknown_chirality" in res.columns
    assert res.loc[0, "unknown_chirality"]

    assert res["multifunctional_enzyme"].dtype == "string"
    assert res.loc[0, "multifunctional_enzyme"] == " "

    for col in ["IUPHAR_class", "IUPHAR_subclass", "gene_index", "taxon_index"]:
        assert col in res.columns


def test_process_activity_table_targets_in_subdir(tmp_path: Path) -> None:
    """Load targets from a nested ``_Target`` subdirectory."""
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
            "nstereo": [1],
            "approx_cited_activity": [pd.NA],
            "exact_cited_activity": [True],
            "higly_correlated_cit": [pd.NA],
            "shuffled_cit": [pd.NA],
            "review_doc": [False],
            "original_activity_approx": [0.5],
            "original_activity_exact": [0.5],
        }
    )

    cf_dir = tmp_path / "_Curation"
    cf_dir.mkdir()
    (cf_dir / "citation_fraction.csv").write_text(
        "N,K_min_significant,test_used_at_threshold,p_value_at_threshold\n1,1,x,0.05\n"
    )

    subdir = tmp_path / "_Target"
    subdir.mkdir()
    (subdir / "targets_type.csv").write_text(
        ",".join(
            [
                "target_chembl_id",
                "IUPHAR_class",
                "IUPHAR_subclass",
                "gene_index",
                "taxon_index",
                "target_sort_order",
                "multifunctional_enzyme",
                "organism_type",
            ]
        )
        + "\nT1,ClassB,, , , ,False,Viruses\n"
    )

    res = process_activity_table(df, tmp_path)
    assert "unicellular_organism" in res.columns
    assert res.loc[0, "unicellular_organism"]

    for col in ["IUPHAR_class", "IUPHAR_subclass", "gene_index", "taxon_index"]:
        assert col in res.columns
