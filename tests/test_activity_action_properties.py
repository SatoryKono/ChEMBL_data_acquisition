from __future__ import annotations

import json

import pandas as pd
import pytest

from library.activity_action_properties import (
    annotate_action_properties,
    build_activity_properties,
    infer_action_type,
)


def _activity_record(**overrides: object) -> dict[str, object]:
    base = {
        "activity_comment": "",
        "assay_description": "",
        "standard_type": "IC50",
        "standard_value": 10.0,
        "standard_relation": "=",
        "standard_units": "nM",
        "relation": "=",
        "units": "nM",
    }
    base.update(overrides)
    return base


def test_infer_action_type_positive_allosteric() -> None:
    features = {"positive": True, "negative": False, "positive_terms": ["pam"]}
    assert infer_action_type(features) == "PAM"


def test_infer_action_type_negative_allosteric() -> None:
    features = {"positive": False, "negative": True, "negative_terms": ["nam"]}
    assert infer_action_type(features) == "NAM"


def test_infer_action_type_conflict_logged(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[dict[str, object]] = []

    def capture(*args, **kwargs) -> None:  # type: ignore[no-untyped-def]
        calls.append({"args": args, "kwargs": kwargs})

    monkeypatch.setattr("library.activity_action_properties.logger.warning", capture)

    result = infer_action_type(
        {
            "positive": True,
            "negative": True,
            "positive_terms": ["pam"],
            "negative_terms": ["nam"],
        }
    )
    assert result is None
    assert calls and calls[0]["args"][0] == "effect_conflict"


def test_build_activity_properties_logs_unmapped_metric(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[dict[str, object]] = []

    def capture(*args, **kwargs) -> None:  # type: ignore[no-untyped-def]
        calls.append({"args": args, "kwargs": kwargs})

    monkeypatch.setattr("library.activity_action_properties.logger.info", capture)

    payload = build_activity_properties(
        _activity_record(standard_type="TGI50", activity_comment="Positive allosteric modulator")
    )
    assert payload["measurement"]["type"] == "TGI50"
    assert calls and calls[0]["args"][0] == "metric_unmapped"


def test_build_activity_properties_filters_empty_values() -> None:
    payload = build_activity_properties(
        _activity_record(
            activity_comment=" ",
            data_validity_comment="",
            assay_description="Calcium flux",
            assay_variant_accession=None,
            assay_variant_mutation="",
            standard_units="",
            qudt_units=None,
        )
    )
    assert "comments" not in payload
    assert payload["assay"] == {"description": "Calcium flux"}
    assert "qudt_units" not in payload["measurement"]


def test_annotate_action_properties_generates_json() -> None:
    df = pd.DataFrame(
        [
            {
                "activity_id": "A1",
                "molecule_chembl_id": "CHEMBL1",
                "assay_chembl_id": "AS1",
                "standard_type": "IC50",
                "standard_value": 1.0,
                "activity_comment": "positive allosteric modulator",
            }
        ]
    )
    result = annotate_action_properties(df)
    assert "activity_properties" in result
    payload = json.loads(result.loc[0, "activity_properties"])
    assert result.loc[0, "action_type"] == "PAM"
    assert payload["measurement"]["value"] == 1.0

