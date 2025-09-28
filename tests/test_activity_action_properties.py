import json

import pandas as pd
import pytest

from library.activity_action_properties import annotate_action_properties
from library.config import ActivityActionTypeCfg, ActivityPropertiesCfg
from scripts.get_activity_data import (
    _extract_effect_features,
    _normalise_mapping,
    apply_activity_annotations,
    build_activity_properties,
    infer_action_type,
)


def _record(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "assay_chembl_id": "AS1",
        "assay_description": "",
        "activity_comment": "",
        "standard_type": "IC50",
        "standard_value": 10.0,
        "standard_relation": "=",
        "standard_units": "nM",
        "relation": "=",
        "units": "nM",
    }
    base.update(overrides)
    return base


def _tracker() -> dict[str, set[str]]:
    return {
        "metric_unmapped": set(),
        "functionality_unmapped": set(),
        "mechanism_unmapped": set(),
        "unexpected_outputs": set(),
        "triage_unmapped": set(),
    }


def test_infer_action_type_metric_precedence() -> None:
    cfg = ActivityActionTypeCfg()
    record = _record(activity_comment="positive allosteric modulator")
    features = _extract_effect_features(record)
    allowlist = {token for value in cfg.allowlist if (token := value.strip().lower())}
    result = infer_action_type(
        record,
        cfg,
        features=features,
        metrics_map=_normalise_mapping(cfg.metrics),
        functionality_map=_normalise_mapping(cfg.functionality),
        mechanism_map=_normalise_mapping(cfg.mechanism),
        allowlist=allowlist,
        trackers=_tracker(),
        positive_label=cfg.positive_label,
        negative_label=cfg.negative_label,
        fallback_label=cfg.fallback,
    )
    assert result == "inhibition"


def test_infer_action_type_allosteric_fallback() -> None:
    cfg = ActivityActionTypeCfg()
    record = _record(
        standard_type="TGI50", activity_comment="positive allosteric modulator"
    )
    features = _extract_effect_features(record)
    allowlist = {token for value in cfg.allowlist if (token := value.strip().lower())}
    trackers = _tracker()
    result = infer_action_type(
        record,
        cfg,
        features=features,
        metrics_map=_normalise_mapping(cfg.metrics),
        functionality_map=_normalise_mapping(cfg.functionality),
        mechanism_map=_normalise_mapping(cfg.mechanism),
        allowlist=allowlist,
        trackers=trackers,
        positive_label=cfg.positive_label,
        negative_label=cfg.negative_label,
        fallback_label=cfg.fallback,
    )
    assert result == cfg.positive_label
    assert "tgi50" in trackers["metric_unmapped"]


def test_build_activity_properties_generates_hash() -> None:
    cfg = ActivityPropertiesCfg()
    record = _record(
        activity_comment="Active",
        data_validity_comment="Triaged",
        assay_description="Calcium flux",
        assay_variant_accession=None,
        assay_variant_mutation="",
    )
    json_text, hash_value = build_activity_properties(
        record,
        cfg,
        features=_extract_effect_features(record),
        metrics_map=_normalise_mapping(ActivityActionTypeCfg().metrics),
        triage_map=_normalise_mapping({"triaged": "triaged"}),
        triage_fields=["data_validity_comment"],
        functionality_fields=["functional_activity"],
        mechanism_fields=["mechanism_of_action"],
        triage_unmapped=set(),
    )
    assert json_text is not None
    payload = json.loads(json_text)
    assert payload["measurement"]["type"] == "IC50"
    assert payload["assay"] == {"assay_chembl_id": "AS1", "description": "Calcium flux"}
    assert hash_value is not None and len(hash_value) == 64


def test_apply_activity_annotations_adds_columns() -> None:
    action_cfg = ActivityActionTypeCfg()
    properties_cfg = ActivityPropertiesCfg()
    frame = pd.DataFrame(
        [
            _record(
                activity_comment="positive allosteric modulator",
                data_validity_comment="triaged",
                standard_type="TGI50",
                standard_value=1.5,
            )
        ]
    )
    annotated = apply_activity_annotations(frame, action_cfg, properties_cfg)
    assert action_cfg.column in annotated
    assert properties_cfg.column in annotated
    assert properties_cfg.hash_column in annotated
    payload = json.loads(annotated.loc[0, properties_cfg.column])
    assert annotated.loc[0, action_cfg.column] == "PAM"
    assert pytest.approx(payload["measurement"]["value"], rel=1e-6) == 1.5


def test_annotate_action_properties_streams_rows(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_to_dict = pd.DataFrame.to_dict
    original_itertuples = pd.DataFrame.itertuples

    def fail_on_records(
        self: pd.DataFrame, orient: str = "dict", *args: object, **kwargs: object
    ):
        if orient == "records":
            raise AssertionError("to_dict with orient='records' should not be used")
        return original_to_dict(self, orient, *args, **kwargs)

    def tracking_itertuples(
        self: pd.DataFrame, index: bool = False, name: str | None = None
    ):
        call_count[0] += 1
        iterator = original_itertuples(self, index=index, name=name)

        def generator():
            for row in iterator:
                yielded_rows[0] += 1
                yield row

        return generator()

    call_count = [0]
    yielded_rows = [0]
    monkeypatch.setattr(pd.DataFrame, "to_dict", fail_on_records)
    monkeypatch.setattr(pd.DataFrame, "itertuples", tracking_itertuples)

    large_frame = pd.DataFrame(
        [
            _record(
                activity_comment="positive allosteric modulator",
                standard_value=index,
            )
            for index in range(1000)
        ]
    )

    annotated = annotate_action_properties(large_frame)

    assert call_count[0] == 1
    assert yielded_rows[0] == len(large_frame)
    assert annotated.shape[0] == len(large_frame)
    assert "activity_properties" in annotated
    assert "action_type" in annotated
