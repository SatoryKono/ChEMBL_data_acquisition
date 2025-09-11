"""Tests for status processing pipeline."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.input_initialisation_library import (
    aggregate_activity,
    build_status_helpers,
    initialize_activity_status,
    initialize_pairs,
    load_status_table,
)


def test_status_pipeline(tmp_path) -> None:
    """Smoke test the status workflow end-to-end."""
    (tmp_path / "status.csv").write_text(
        "\n".join(
            [
                "status,condition_field,condition_value,order,score",
                "good,null,null,0,1",
                "warning,review,True,1,0",
                "bad,high_citation_rate,True,2,-1",
            ]
        )
    )
    status_df = load_status_table(tmp_path)
    api = build_status_helpers(status_df)
    assert api.status_list == ["good", "warning", "bad"]
    assert api.conditions == ["review", "high_citation_rate"]
    assert api.pair("warning", "bad") == "warning"
    assert api.Next("bad") == "bad"

    activity = pd.DataFrame(
        {
            "activity_id": [1, 2],
            "assay_id": ["A1", "A1"],
            "document_id": ["D1", "D1"],
            "testitem_id": ["T1", "T2"],
            "target_id": ["TG1", "TG1"],
            "mesurement_type": ["IC50", "IC50"],
            "high_citation_rate": [False, False],
            "unicellular_organism": [False, False],
            "review": [False, True],
            "rounded_data_citation": [False, False],
            "shuffled_assay": [False, False],
            "higly_correlated_assay": [False, False],
            "exact_data_citation": [False, False],
            "multmol_assay": [False, False],
            "multifunctional_enzyme": [False, False],
            "unknown_chirality": [False, False],
        }
    )
    activity = initialize_activity_status(activity, api)
    assert list(activity["Filtered.init"]) == ["bad", "warning"]

    pair_df = pd.DataFrame(
        {
            "activity_chembl_id1": [1],
            "activity_chembl_id2": [2],
            "testitem_id": ["T1"],
            "target_id": ["TG1"],
            "mesurement_type": ["IC50"],
            "independent_IC50": [1],
            "non_independent_IC50": [0],
            "independent_Ki": [0],
            "non_independent_Ki": [0],
        }
    )
    pair_df = initialize_pairs(pair_df, activity, api)
    assert pair_df.loc[0, "Filtered"] == "warning"
    agg = aggregate_activity(pair_df, activity, api)
    assert agg["activity"].set_index("activity_id").loc[1, "Filtered.new"] == "bad"
    assert agg["assay"].loc[0, "independent_IC50"] == 1


def test_load_status_table_skips_empty_rows(tmp_path: Path) -> None:
    """``load_status_table`` should ignore completely empty CSV lines."""
    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n"
        "good,null,null,0,1\n"
        "\n"
        "bad,null,null,1,2\n"
    )
    df = load_status_table(tmp_path)
    assert df["status"].tolist() == ["good", "bad"]


def test_aggregate_activity_handles_missing_metrics(tmp_path: Path) -> None:
    """``aggregate_activity`` should default missing metric columns to zeros."""
    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n" "good,null,null,0,0\n"
    )
    status_df = load_status_table(tmp_path)
    api = build_status_helpers(status_df)
    activity = pd.DataFrame(
        {
            "activity_id": [1, 2],
            "assay_id": ["A1", "A1"],
            "document_id": ["D1", "D1"],
            "testitem_id": ["T1", "T2"],
            "target_id": ["TG1", "TG1"],
            "mesurement_type": ["IC50", "IC50"],
            "high_citation_rate": [False, False],
            "unicellular_organism": [False, False],
        }
    )
    activity = initialize_activity_status(activity, api)
    pair_df = pd.DataFrame(
        {
            "activity_chembl_id1": [1],
            "activity_chembl_id2": [2],
            "testitem_id": ["T1"],
            "target_id": ["TG1"],
            "mesurement_type": ["IC50"],
        }
    )
    pair_df = initialize_pairs(pair_df, activity, api)
    agg = aggregate_activity(pair_df, activity, api)
    activity_status = agg["activity"].set_index("activity_id")
    for col in [
        "independent_IC50",
        "non_independent_IC50",
        "independent_Ki",
        "non_independent_Ki",
    ]:
        assert col in activity_status.columns
        assert activity_status[col].eq(0).all()


def test_aggregate_activity_handles_missing_testitem_id(tmp_path: Path) -> None:
    """Missing ``testitem_id`` should not break aggregation."""
    (tmp_path / "status.csv").write_text(
        "status,condition_field,condition_value,order,score\n" "good,null,null,0,0\n"
    )
    status_df = load_status_table(tmp_path)
    api = build_status_helpers(status_df)
    # activity table deliberately lacks ``testitem_id``
    activity = pd.DataFrame(
        {
            "activity_id": [1, 2],
            "assay_id": ["A1", "A1"],
            "document_id": ["D1", "D1"],
            "target_id": ["TG1", "TG1"],
            "mesurement_type": ["IC50", "IC50"],
            "high_citation_rate": [False, False],
            "unicellular_organism": [False, False],
        }
    )
    activity = initialize_activity_status(activity, api)
    pair_df = pd.DataFrame(
        {
            "activity_chembl_id1": [1],
            "activity_chembl_id2": [2],
            "independent_IC50": [1],
            "non_independent_IC50": [0],
            "independent_Ki": [0],
            "non_independent_Ki": [0],
        }
    )
    pair_df = initialize_pairs(pair_df, activity, api)
    agg = aggregate_activity(pair_df, activity, api)
    assert agg["system"].empty
    assert agg["testitem"].empty
    assert agg["target"].set_index("target_id").loc["TG1", "independent_IC50"] == 1
