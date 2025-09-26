"""Derive action annotations and structured metadata for activity records."""

from __future__ import annotations

from collections.abc import Mapping
import json
from typing import Any

import pandas as pd

from .log import logger

# Keywords describing allosteric modulation effects.  The lists intentionally
# focus on clear-cut textual cues that commonly appear in ChEMBL annotations.
_POSITIVE_KEYWORDS = {
    "positive allosteric",  # explicit mention in comments/assay descriptions
    "pam",  # abbreviation frequently used by submitters
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

# Mapping between measurement types and a coarse category used for downstream
# aggregation.  The mapping is intentionally small; values falling outside the
# dictionary keep their ``standard_type`` but trigger a log entry so that data
# engineers notice new metrics.
_MEASUREMENT_KIND = {
    "IC50": "inhibition",
    "EC50": "activation",
    "AC50": "activation",
    "Ki": "binding",
    "Kd": "binding",
}


def _has_value(value: Any) -> bool:
    """Return ``True`` if *value* should be preserved in the JSON payload."""

    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, set)):
        return any(_has_value(item) for item in value)
    if isinstance(value, Mapping):
        return any(_has_value(item) for item in value.values())
    try:
        return not pd.isna(value)
    except TypeError:  # pragma: no cover - non-numeric objects
        return True


def _clean_structure(data: Any) -> Any:
    """Recursively remove empty values from nested ``dict``/``list`` objects."""

    if isinstance(data, Mapping):
        cleaned = {
            key: _clean_structure(value)
            for key, value in data.items()
            if _has_value(value)
        }
        return cleaned if cleaned else None
    if isinstance(data, list):
        cleaned_list = [item for item in (_clean_structure(v) for v in data) if _has_value(item)]
        return cleaned_list if cleaned_list else None
    if isinstance(data, tuple):
        cleaned_tuple = tuple(
            item for item in (_clean_structure(v) for v in data) if _has_value(item)
        )
        return cleaned_tuple if cleaned_tuple else None
    return data


def _coalesce(*values: Any) -> Any:
    """Return the first non-empty entry from ``values``."""

    for value in values:
        if _has_value(value):
            return value
    return None


def _collect_text(record: Mapping[str, Any]) -> str:
    """Concatenate textual annotations for keyword detection."""

    parts: list[str] = []
    for key in (
        "activity_comment",
        "data_validity_comment",
        "data_validity_description",
        "assay_description",
    ):
        value = record.get(key)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                parts.append(stripped)
    return " ".join(parts).lower()


def _extract_effect_features(record: Mapping[str, Any]) -> dict[str, Any]:
    """Return effect flags inferred from textual annotations in ``record``."""

    text = _collect_text(record)
    positive_hits = sorted({kw for kw in _POSITIVE_KEYWORDS if kw in text})
    negative_hits = sorted({kw for kw in _NEGATIVE_KEYWORDS if kw in text})
    is_allosteric = bool(positive_hits or negative_hits or any(kw in text for kw in _ALLOSTERIC_KEYWORDS))
    return {
        "allosteric": is_allosteric,
        "positive": bool(positive_hits),
        "negative": bool(negative_hits),
        "positive_terms": positive_hits,
        "negative_terms": negative_hits,
    }


def infer_action_type(features: Mapping[str, Any]) -> str | None:
    """Infer the high-level action type from *features* flags.

    ``features`` is expected to contain boolean flags ``positive`` and
    ``negative``.  When both are present the effect is ambiguous and the
    function returns ``None`` after logging a warning.  Positive hits map to the
    ``PAM`` label and negative hits map to ``NAM``.
    """

    positive = bool(features.get("positive"))
    negative = bool(features.get("negative"))
    if positive and negative:
        logger.warning(
            "effect_conflict",
            positive_terms=features.get("positive_terms", []),
            negative_terms=features.get("negative_terms", []),
        )
        return None
    if positive:
        return "PAM"
    if negative:
        return "NAM"
    return None


def build_activity_properties(record: Mapping[str, Any]) -> dict[str, Any]:
    """Return a structured metadata payload for ``record`` suitable for JSON."""

    features = _extract_effect_features(record)
    measurement_type = _coalesce(record.get("standard_type"), record.get("type"))
    units = _coalesce(record.get("standard_units"), record.get("units"))
    relation = _coalesce(record.get("standard_relation"), record.get("relation"))
    measurement: dict[str, Any] = {
        "type": measurement_type,
        "value": record.get("standard_value"),
        "relation": relation,
        "units": units,
        "lower_value": record.get("lower_value"),
        "upper_value": record.get("upper_value"),
        "pchembl_value": record.get("pchembl_value"),
        "qudt_units": record.get("qudt_units"),
    }
    kind = _MEASUREMENT_KIND.get(measurement_type)
    if measurement_type and kind is None:
        logger.info("metric_unmapped", standard_type=measurement_type)
    if kind:
        measurement["kind"] = kind
    measurement = _clean_structure(measurement) or {}

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
        value.strip()
        for value in (
            record.get("activity_comment"),
            record.get("data_validity_comment"),
            record.get("data_validity_description"),
        )
        if isinstance(value, str) and value.strip()
    ]

    payload: dict[str, Any] = {"measurement": measurement}
    if assay_info:
        payload["assay"] = assay_info
    if comments:
        payload["comments"] = comments
    if features:
        payload["effect_features"] = features
    return payload


def annotate_action_properties(df: pd.DataFrame) -> pd.DataFrame:
    """Attach ``activity_properties`` JSON and ``action_type`` to ``df``."""

    if df.empty:
        result = df.copy()
        result["activity_properties"] = pd.Series(dtype=object)
        result["action_type"] = pd.Series(dtype=object)
        return result

    properties: list[str] = []
    action_types: list[str | None] = []
    for record in df.to_dict(orient="records"):
        payload = build_activity_properties(record)
        features = payload.get("effect_features", {})
        action_types.append(infer_action_type(features))
        properties.append(json.dumps(payload, sort_keys=True))

    result = df.copy()
    result["activity_properties"] = properties
    result["action_type"] = action_types
    return result


__all__ = [
    "annotate_action_properties",
    "build_activity_properties",
    "infer_action_type",
]

