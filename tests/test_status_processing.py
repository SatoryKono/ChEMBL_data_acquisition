"""Tests for status processing pipeline."""

from __future__ import annotations

import pandas as pd

from library.input_initialisation_library import (
    aggregate_activity,
    build_status_helpers,
    initialize_activity_status,
    initialize_pairs,
    load_status_table,
)


def test_status_pipeline(tmp_path):
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
            "activity_id1": [1],
            "activity_id2": [2],
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
