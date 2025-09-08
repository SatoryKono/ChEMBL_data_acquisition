from __future__ import annotations

import pandas as pd
import sys
from pathlib import Path

import pytest

from library.input_initialisation_library import (
    _ensure_openpyxl,
    EntityName,
    append_entities,
    build_combined_tables,
    unify_dtypes,
    save_tables,
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


def test_build_combined_tables_drops_activity_cols() -> None:
    same: dict[EntityName, pd.DataFrame] = {
        "activity": pd.DataFrame({"id": [1], "Column1": ["a"]}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
    }
    all_: dict[EntityName, pd.DataFrame] = {
        "activity": pd.DataFrame({"id": [2], "Column1": ["b"]}),
        "assay": pd.DataFrame(),
        "document": pd.DataFrame(),
        "target": pd.DataFrame(),
        "testitem": pd.DataFrame(),
    }
    combined = build_combined_tables(same, all_)
    assert "Column1" not in combined["activity"].columns
    assert len(combined["activity"]) == 2


def test_save_tables_writes_files(tmp_path: Path) -> None:
    tables: dict[EntityName, pd.DataFrame] = {
        "activity": pd.DataFrame({"id": [1]}),
        "assay": pd.DataFrame({"id": [2]}),
        "document": pd.DataFrame({"id": [3]}),
        "target": pd.DataFrame({"id": [4]}),
        "testitem": pd.DataFrame({"id": [5]}),
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
