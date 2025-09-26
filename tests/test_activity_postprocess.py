"""Tests for activity post-processing helpers."""

from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest

from library.config import ActivityActionCfg, ActivityPropertiesCfg
from scripts.get_activity_data import build_activity_properties, infer_action_type


def test_infer_action_type_priority_and_logging(monkeypatch: pytest.MonkeyPatch) -> None:
    """Metric mapping takes precedence and unmapped values are logged once."""

    cfg = ActivityActionCfg(
        metrics={"IC50": "inhibition"},
        allosteric={"positive allosteric modulator": "modulation_positive"},
        functionality={"agonist activity": "agonist"},
        mechanism={"inhibitor": "inhibition"},
        fallback="unknown",
    )
    df = pd.DataFrame(
        [
            {
                "standard_type": "IC50",
                "activity_comment": "Positive Allosteric Modulator",
                "bao_label": None,
            },
            {
                "standard_type": None,
                "bao_label": "Agonist Activity",
            },
            {
                "standard_type": None,
                "activity_comment": "Inhibitor",
            },
            {
                "standard_type": None,
                "activity_comment": "Unknown Effect",
            },
        ]
    )

    captured: list[dict[str, object]] = []

    def fake_warning(event: str, **extra: object) -> None:
        if event == "action_type_unmapped":
            captured.append(extra)

    monkeypatch.setattr(
        "scripts.get_activity_data.logger.warning",
        fake_warning,
    )
    result = infer_action_type(df, cfg)

    assert result.tolist() == [
        "inhibition",
        "agonist",
        "inhibition",
        "unknown",
    ]
    categories = {entry.get("category") for entry in captured}
    assert {"allosteric", "mechanism"}.issubset(categories)
    assert any("unknown effect" in entry.get("values", []) for entry in captured)


def test_build_activity_properties_serialization(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Properties JSON honours allowlist, mapping and hashing rules."""

    cfg = ActivityPropertiesCfg(
        allowlist=["activity_comment", "data_validity_description"],
        triage={"data_validity_description": {"manually validated": "validated"}},
        hash_fields=["activity_comment", "data_validity_description"],
    )
    df = pd.DataFrame(
        [
            {
                "activity_comment": " Strong effect ",
                "data_validity_description": "Manually Validated",
            },
            {
                "activity_comment": "   ",
                "data_validity_description": "Unknown",
            },
            {
                "activity_comment": None,
                "data_validity_description": None,
            },
        ]
    )

    triage_events: list[dict[str, object]] = []

    def fake_warning(event: str, **extra: object) -> None:
        if event == "activity_properties_unmapped_triage":
            triage_events.append(extra)

    monkeypatch.setattr(
        "scripts.get_activity_data.logger.warning",
        fake_warning,
    )
    props, hashes = build_activity_properties(df, cfg)

    assert props.tolist() == [
        "{\"activity_comment\":\"Strong effect\",\"data_validity_description\":\"validated\"}",
        "{\"data_validity_description\":\"Unknown\"}",
        None,
    ]

    expected_first = hashlib.sha256(
        json.dumps(
            {
                "activity_comment": "Strong effect",
                "data_validity_description": "validated",
            },
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    expected_second = hashlib.sha256(
        json.dumps(
            {"data_validity_description": "Unknown"},
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()
    assert hashes is not None
    assert hashes.tolist() == [expected_first, expected_second, None]

    assert triage_events
    assert triage_events[0].get("column") == "data_validity_description"
    assert "unknown" in triage_events[0].get("values", [])
