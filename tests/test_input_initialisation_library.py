from __future__ import annotations

import pandas as pd

from library.input_initialisation_library import append_entities, unify_dtypes


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
