"""Helpers for deriving activity annotations and metadata."""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from hashlib import sha256
from typing import Any

import pandas as pd

from ...common.log import logger
from ...config import ActivityActionTypeCfg, ActivityPropertiesCfg

__all__ = [
    "apply_activity_annotations",
    "build_activity_properties",
    "infer_action_type",
    "normalise_mapping",
    "extract_effect_features",
]

_POSITIVE_KEYWORDS = {
    "positive allosteric",
    "pam",
}
_NEGATIVE_KEYWORDS = {
    "negative allosteric",
    "nam",
}
_ALLOSTERIC_KEYWORDS = {
    "allosteric modulator",
    "allosteric inhibitor",
    "allosteric activator",
}
_METRIC_FIELDS: tuple[str, ...] = ("standard_type", "type")
_COMMENT_FIELDS: tuple[str, ...] = (
    "activity_comment",
    "data_validity_comment",
    "data_validity_description",
    "assay_description",
)


def _normalise_text(value: object) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    text = str(value).strip()
    return text or None


def _normalise_token(value: object) -> str | None:
    if value is None:
        return None
    if isinstance(value, str):
        trimmed = value.strip()
        return trimmed.lower() or None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    text = str(value).strip()
    return text.lower() or None


def _normalise_output(value: object) -> str | None:
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    text = str(value).strip()
    return text or None


def _has_value(value: Any) -> bool:
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, list | tuple | set):
        return any(_has_value(item) for item in value)
    if isinstance(value, Mapping):
        return any(_has_value(item) for item in value.values())
    try:
        return not pd.isna(value)
    except TypeError:
        return True


def _clean_structure(data: Any) -> Any:
    if isinstance(data, Mapping):
        cleaned = {
            key: _clean_structure(value)
            for key, value in data.items()
            if _has_value(value)
        }
        return cleaned or None
    if isinstance(data, list):
        cleaned_list = [
            item
            for item in (_clean_structure(value) for value in data)
            if _has_value(item)
        ]
        return cleaned_list or None
    if isinstance(data, tuple):
        cleaned_tuple = tuple(
            item
            for item in (_clean_structure(value) for value in data)
            if _has_value(item)
        )
        return cleaned_tuple or None
    return data


def _collect_text(record: Mapping[str, Any]) -> str:
    parts: list[str] = []
    for field in _COMMENT_FIELDS:
        value = record.get(field)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                parts.append(stripped)
    return " ".join(parts).lower()


def extract_effect_features(record: Mapping[str, Any]) -> dict[str, Any]:
    text = _collect_text(record)
    positive_hits = sorted({kw for kw in _POSITIVE_KEYWORDS if kw in text})
    negative_hits = sorted({kw for kw in _NEGATIVE_KEYWORDS if kw in text})
    is_allosteric = bool(
        positive_hits
        or negative_hits
        or any(keyword in text for keyword in _ALLOSTERIC_KEYWORDS)
    )
    return {
        "allosteric": is_allosteric,
        "positive": bool(positive_hits),
        "negative": bool(negative_hits),
        "positive_terms": positive_hits,
        "negative_terms": negative_hits,
    }


def _first_non_empty(values: Iterable[object]) -> str | None:
    for value in values:
        text = _normalise_output(value)
        if text:
            return text
    return None


def _lookup_from_fields(
    record: Mapping[str, Any],
    fields: Sequence[str],
    mapping: Mapping[str, str],
    *,
    unmapped_inputs: set[str],
    unexpected_outputs: set[str],
    allowlist: set[str],
) -> str | None:
    for field in fields:
        token = _normalise_token(record.get(field))
        if token is None:
            continue
        mapped = mapping.get(token)
        if mapped is None:
            if token not in allowlist:
                unmapped_inputs.add(token)
            continue
        cleaned = _normalise_output(mapped)
        if not cleaned:
            continue
        if allowlist and cleaned.lower() not in allowlist:
            unexpected_outputs.add(cleaned)
            continue
        return cleaned
    return None


def infer_action_type(
    record: Mapping[str, Any],
    cfg: ActivityActionTypeCfg,
    *,
    features: Mapping[str, Any],
    metrics_map: Mapping[str, str],
    functionality_map: Mapping[str, str],
    mechanism_map: Mapping[str, str],
    allowlist: set[str],
    trackers: dict[str, set[str]],
    positive_label: str | None,
    negative_label: str | None,
    fallback_label: str | None,
) -> str | None:
    metric_value = _lookup_from_fields(
        record,
        _METRIC_FIELDS,
        metrics_map,
        unmapped_inputs=trackers["metric_unmapped"],
        unexpected_outputs=trackers["unexpected_outputs"],
        allowlist=allowlist,
    )
    if metric_value:
        return metric_value

    if features.get("positive") and positive_label:
        if allowlist and positive_label.lower() not in allowlist:
            trackers["unexpected_outputs"].add(positive_label)
        else:
            return positive_label

    if features.get("negative") and negative_label:
        if allowlist and negative_label.lower() not in allowlist:
            trackers["unexpected_outputs"].add(negative_label)
        else:
            return negative_label

    functionality_value = _lookup_from_fields(
        record,
        cfg.functionality_fields,
        functionality_map,
        unmapped_inputs=trackers["functionality_unmapped"],
        unexpected_outputs=trackers["unexpected_outputs"],
        allowlist=allowlist,
    )
    if functionality_value:
        return functionality_value

    mechanism_value = _lookup_from_fields(
        record,
        cfg.mechanism_fields,
        mechanism_map,
        unmapped_inputs=trackers["mechanism_unmapped"],
        unexpected_outputs=trackers["unexpected_outputs"],
        allowlist=allowlist,
    )
    if mechanism_value:
        return mechanism_value

    if fallback_label:
        if allowlist and fallback_label.lower() not in allowlist:
            trackers["unexpected_outputs"].add(fallback_label)
        else:
            return fallback_label
    return None


def build_activity_properties(
    record: Mapping[str, Any],
    properties_cfg: ActivityPropertiesCfg,
    *,
    features: Mapping[str, Any],
    metrics_map: Mapping[str, str],
    triage_map: Mapping[str, str],
    triage_fields: Sequence[str],
    functionality_fields: Sequence[str],
    mechanism_fields: Sequence[str],
    triage_unmapped: set[str],
) -> tuple[str | None, str | None]:
    measurement_type = _first_non_empty(record.get(field) for field in _METRIC_FIELDS)
    measurement_kind = None
    if measurement_type:
        token = _normalise_token(measurement_type)
        if token:
            measurement_kind = _normalise_output(metrics_map.get(token))
    measurement = _clean_structure(
        {
            "type": measurement_type,
            "kind": measurement_kind,
            "value": record.get("standard_value"),
            "relation": _first_non_empty(
                [record.get("standard_relation"), record.get("relation")]
            ),
            "units": _first_non_empty(
                [record.get("standard_units"), record.get("units")]
            ),
            "lower_value": record.get("lower_value"),
            "upper_value": record.get("upper_value"),
            "pchembl_value": record.get("pchembl_value"),
            "qudt_units": record.get("qudt_units"),
        }
    )

    assay_info = _clean_structure(
        {
            "assay_chembl_id": record.get("assay_chembl_id"),
            "description": record.get("assay_description"),
            "bao_label": record.get("bao_label"),
            "bao_format": record.get("bao_format"),
            "variant": {
                "accession": record.get("assay_variant_accession"),
                "mutation": record.get("assay_variant_mutation"),
            },
        }
    )

    comments = [
        value
        for value in (
            _normalise_output(record.get("activity_comment")),
            _normalise_output(record.get("data_validity_comment")),
            _normalise_output(record.get("data_validity_description")),
        )
        if value
    ]

    triage_payload: dict[str, str] | None = None
    for field in triage_fields:
        raw_value = record.get(field)
        token = _normalise_token(raw_value)
        if token is None:
            continue
        mapped = triage_map.get(token)
        if mapped is None:
            triage_unmapped.add(token)
            continue
        label = _normalise_output(mapped)
        source_value = _normalise_output(raw_value)
        if label is None or source_value is None:
            continue
        triage_payload = {
            "label": label,
            "source_field": field,
            "source_value": source_value,
        }
        break

    functionality_value = _first_non_empty(
        record.get(field) for field in functionality_fields
    )
    mechanism_value = _first_non_empty(record.get(field) for field in mechanism_fields)

    payload = {
        "measurement": measurement,
        "assay": assay_info,
        "comments": comments or None,
        "effect_features": features or None,
        "triage": triage_payload,
        "functionality": functionality_value,
        "mechanism": mechanism_value,
    }

    allowed_keys = {key for key in properties_cfg.allowlist}
    filtered = {key: value for key, value in payload.items() if key in allowed_keys}
    cleaned = _clean_structure(filtered)
    if not cleaned:
        return None, None

    json_text = json.dumps(
        cleaned, ensure_ascii=False, separators=(",", ":"), sort_keys=True
    )
    hash_value = None
    if properties_cfg.hash_column:
        hash_value = sha256(json_text.encode("utf-8")).hexdigest()
    return json_text, hash_value


def normalise_mapping(mapping: Mapping[str, str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, value in mapping.items():
        token = _normalise_token(key)
        cleaned_value = _normalise_output(value)
        if token and cleaned_value:
            result[token] = cleaned_value
    return result


def apply_activity_annotations(
    df: pd.DataFrame,
    action_cfg: ActivityActionTypeCfg,
    properties_cfg: ActivityPropertiesCfg,
) -> pd.DataFrame:
    if df.empty:
        result = df.copy()
        if action_cfg.enabled:
            result[action_cfg.column] = pd.Series(dtype="string")
        if properties_cfg.enabled:
            result[properties_cfg.column] = pd.Series(dtype="string")
        if properties_cfg.hash_column:
            result[properties_cfg.hash_column] = pd.Series(dtype="string")
        return result

    metrics_map = normalise_mapping(action_cfg.metrics)
    functionality_map = normalise_mapping(action_cfg.functionality)
    mechanism_map = normalise_mapping(action_cfg.mechanism)
    triage_map = normalise_mapping(action_cfg.triages)
    allowlist = {
        token
        for value in action_cfg.allowlist
        if (token := _normalise_token(value)) is not None
    }

    positive_label = _normalise_output(action_cfg.positive_label)
    negative_label = _normalise_output(action_cfg.negative_label)
    fallback_label = _normalise_output(action_cfg.fallback)

    trackers: dict[str, set[str]] = {
        "metric_unmapped": set(),
        "functionality_unmapped": set(),
        "mechanism_unmapped": set(),
        "unexpected_outputs": set(),
        "triage_unmapped": set(),
    }

    action_values: list[str | None] = []
    property_values: list[str | None] = []
    hash_values: list[str | None] = []

    records = df.to_dict(orient="records")
    for record in records:
        features = extract_effect_features(record)
        action = infer_action_type(
            record,
            action_cfg,
            features=features,
            metrics_map=metrics_map,
            functionality_map=functionality_map,
            mechanism_map=mechanism_map,
            allowlist=allowlist,
            trackers=trackers,
            positive_label=positive_label,
            negative_label=negative_label,
            fallback_label=fallback_label,
        )
        action_values.append(action)

        if properties_cfg.enabled or properties_cfg.hash_column:
            json_value, hash_value = build_activity_properties(
                record,
                properties_cfg,
                features=features,
                metrics_map=metrics_map,
                triage_map=triage_map,
                triage_fields=action_cfg.triage_fields,
                functionality_fields=action_cfg.functionality_fields,
                mechanism_fields=action_cfg.mechanism_fields,
                triage_unmapped=trackers["triage_unmapped"],
            )
        else:
            json_value, hash_value = None, None
        property_values.append(json_value)
        hash_values.append(hash_value)

    result = df.copy()
    action_series = pd.Series(action_values, index=result.index, dtype="string")
    if action_cfg.enabled:
        result[action_cfg.column] = action_series
    elif action_cfg.column in result.columns:
        result = result.drop(columns=[action_cfg.column])
    if action_cfg.log_missing:
        missing = int(action_series.isna().sum())
        if missing:
            logger.warning("activity_action_type_missing", rows=missing)
    if action_cfg.log_distribution:
        counts = action_series.dropna().value_counts()
        if not counts.empty:
            logger.info(
                "activity_action_type_distribution",
                counts={str(key): int(value) for key, value in sorted(counts.items())},
            )

    if properties_cfg.enabled:
        properties_series = pd.Series(
            property_values, index=result.index, dtype="string"
        )
        result[properties_cfg.column] = properties_series
        if properties_cfg.log_missing:
            missing = int(properties_series.isna().sum())
            if missing:
                logger.warning("activity_properties_missing", rows=missing)
    elif properties_cfg.drop_source_column and properties_cfg.column in result.columns:
        result = result.drop(columns=[properties_cfg.column])
    if properties_cfg.hash_column:
        hash_series = pd.Series(hash_values, index=result.index, dtype="string")
        result[properties_cfg.hash_column] = hash_series

    if trackers["metric_unmapped"]:
        logger.warning(
            "activity_action_type_unmapped_metric",
            values=sorted(trackers["metric_unmapped"]),
        )
    if trackers["functionality_unmapped"]:
        logger.warning(
            "activity_action_type_unmapped_functionality",
            values=sorted(trackers["functionality_unmapped"]),
        )
    if trackers["mechanism_unmapped"]:
        logger.warning(
            "activity_action_type_unmapped_mechanism",
            values=sorted(trackers["mechanism_unmapped"]),
        )
    if trackers["unexpected_outputs"]:
        logger.warning(
            "activity_action_type_unexpected_output",
            values=sorted(trackers["unexpected_outputs"]),
        )
    if trackers["triage_unmapped"]:
        logger.info(
            "activity_properties_unmapped_triage",
            values=sorted(trackers["triage_unmapped"]),
        )

    return result
