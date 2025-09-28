"""Command line interface for retrieving ChEMBL activity data."""

from __future__ import annotations

import json
import re
import sys
from hashlib import sha256
from pathlib import Path

import argparse
from collections.abc import Iterable, Mapping, Sequence
from itertools import islice
from typing import Any

import numpy as np
import pandas as pd
import requests
from pandera.errors import SchemaErrors

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)

from library.utils import bootstrap

bootstrap.ensure_project_root()

from library import chembl_library as cl
from library import cli
from library import io
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import (
    ActivityActionTypeCfg,
    ActivityBoundsCfg,
    ActivityPropertiesCfg,
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from library.validation import validate_activities
from schemas import ActivitiesSchema, normalize_activities

_UNCERTAINTY_PATTERN = re.compile(
    r"^\s*(?P<value>-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*(?:±|\+/-|\+-)\s*(?P<delta>\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*$"
)

_NONNEGATIVE_TYPE_TOKENS = {
    "ic50",
    "ec50",
    "ac50",
    "gi50",
    "xc50",
    "dc50",
    "kd",
    "ki",
    "potency",
    "mic",
    "activity",
    "lc50",
    "ld50",
    "ed50",
}

_NONNEGATIVE_TYPE_EXCLUSIONS = {
    "delta",
    "dg",
    "enthalpy",
    "gibbs",
    "log",
}

_CONCENTRATION_UNITS = {
    "m",
    "mm",
    "um",
    "nm",
    "pm",
    "fm",
    "mol/l",
    "mol l-1",
    "mmol/l",
    "umol/l",
    "nmol/l",
    "pmol/l",
    "fmol/l",
}

_CONCENTRATION_KEYWORDS = [
    "mol/l",
    "mol l-1",
    "molar",
    "mol per litre",
    "mg/ml",
    "ug/ml",
    "µg/ml",
    "ng/ml",
    "pg/ml",
    "fg/ml",
    "mg/l",
    "ug/l",
    "µg/l",
    "ng/l",
    "pg/l",
    "fg/l",
    "%",
]


def _to_numeric(df: pd.DataFrame, column: str) -> pd.Series:
    """Return *column* as a numeric series with ``NaN`` for invalid entries."""

    if column not in df:
        return pd.Series(np.nan, index=df.index, dtype="float64")
    return pd.to_numeric(df[column], errors="coerce")


def _relation_series(df: pd.DataFrame) -> pd.Series:
    """Return relation values preferring standardised variants."""

    if "standard_relation" in df and "relation" in df:
        return df["standard_relation"].fillna(df["relation"])
    if "standard_relation" in df:
        return df["standard_relation"]
    if "relation" in df:
        return df["relation"]
    return pd.Series([None] * len(df), index=df.index, dtype="object")


def _normalize_relation(value: object) -> str | None:
    """Return a canonical representation of *value* or ``None`` for blanks."""

    if value is None or pd.isna(value):
        return None
    if isinstance(value, str):
        trimmed = value.strip()
        if not trimmed:
            return None
        lowered = trimmed.lower()
        if lowered in {"between", "range"}:
            return lowered
        if trimmed in {"=", "=="}:
            return "="
        if lowered in {"~", "≈"}:
            return "~"
        if trimmed in {">", ">="} or lowered in {">", ">="}:
            return ">="
        if trimmed in {"<", "<="} or lowered in {"<", "<="}:
            return "<="
        return trimmed
    return str(value)


def _needs_nonnegative_clamp(std_type: object, std_units: object) -> bool:
    """Return ``True`` if values should be clamped to the non-negative domain."""

    type_flag = False
    if isinstance(std_type, str):
        lowered = std_type.strip().lower()
        if not any(excl in lowered for excl in _NONNEGATIVE_TYPE_EXCLUSIONS):
            if lowered.startswith("p") and lowered not in {"potency"}:
                type_flag = False
            else:
                type_flag = any(token in lowered for token in _NONNEGATIVE_TYPE_TOKENS)
    unit_flag = False
    if isinstance(std_units, str):
        lowered_units = std_units.strip().lower()
        if lowered_units in _CONCENTRATION_UNITS:
            unit_flag = True
        elif any(keyword in lowered_units for keyword in _CONCENTRATION_KEYWORDS):
            unit_flag = True
    return type_flag or unit_flag


def _apply_uncertainty_bounds(
    df: pd.DataFrame,
    lower: pd.Series,
    upper: pd.Series,
    std_values: pd.Series,
) -> tuple[pd.Series, pd.Series]:
    """Return bounds populated from ``±`` expressions in ``standard_text_value``."""

    text = df.get("standard_text_value")
    if text is None:
        return lower, upper
    matches = text.astype("string").str.extract(_UNCERTAINTY_PATTERN)
    if matches.empty:
        return lower, upper
    base = pd.to_numeric(matches["value"], errors="coerce")
    delta = pd.to_numeric(matches["delta"], errors="coerce")
    mask = lower.isna() & upper.isna() & base.notna() & delta.notna()
    if not mask.any():
        return lower, upper
    lower.loc[mask] = base.loc[mask] - delta.loc[mask]
    upper.loc[mask] = base.loc[mask] + delta.loc[mask]
    mismatch = (
        mask
        & std_values.notna()
        & ~np.isclose(std_values.loc[mask], base.loc[mask], equal_nan=True)
    )
    if mismatch.any():
        logger.warning(
            "activity_bounds_uncertainty_mismatch",
            rows=int(mismatch.sum()),
        )
    return lower, upper


def compute_activity_bounds(
    df: pd.DataFrame, bounds_cfg: ActivityBoundsCfg
) -> pd.DataFrame:
    """Derive ``lower_value`` and ``upper_value`` columns based on configuration."""

    if df.empty:
        df = df.copy()
        df["lower_value"] = pd.Series(dtype="float64")
        df["upper_value"] = pd.Series(dtype="float64")
        return df

    result = df.copy()
    lower = pd.Series(np.nan, index=result.index, dtype="float64")
    upper = pd.Series(np.nan, index=result.index, dtype="float64")

    std_lower = _to_numeric(result, "standard_lower_value")
    std_upper = _to_numeric(result, "standard_upper_value")
    std_values = _to_numeric(result, "standard_value")

    lower = std_lower.combine_first(lower)
    upper = std_upper.combine_first(upper)

    range_mask = std_values.notna() & std_upper.notna()
    if range_mask.any():
        range_lower = pd.Series(np.minimum(std_values, std_upper), index=result.index)
        range_upper = pd.Series(np.maximum(std_values, std_upper), index=result.index)
        lower_fill = range_mask & lower.isna()
        upper_fill = range_mask & upper.isna()
        if lower_fill.any():
            lower.loc[lower_fill] = range_lower.loc[lower_fill]
        if upper_fill.any():
            upper.loc[upper_fill] = range_upper.loc[upper_fill]

    if bounds_cfg.enable_from_uncertainty:
        lower, upper = _apply_uncertainty_bounds(result, lower, upper, std_values)

    relations = _relation_series(result).map(_normalize_relation)
    known_relations = {"=", "~", ">=", "<=", "between", "range"}

    if bounds_cfg.enable_from_relation:
        raw_values = _to_numeric(result, "value")

        required_mask = relations.notna() & std_values.isna()
        conversion_mask = required_mask & raw_values.notna()
        if conversion_mask.any():
            logger.warning(
                "activity_bounds_missing_standard_value",
                rows=int(conversion_mask.sum()),
            )

        eq_mask = relations.isin({"=", "~"}) & std_values.notna()
        if eq_mask.any():
            lower_eq = eq_mask & lower.isna()
            upper_eq = eq_mask & upper.isna()
            lower.loc[lower_eq] = std_values.loc[lower_eq]
            upper.loc[upper_eq] = std_values.loc[upper_eq]

        ge_mask = (relations == ">=") & std_values.notna() & lower.isna()
        if ge_mask.any():
            lower.loc[ge_mask] = std_values.loc[ge_mask]

        le_mask = (relations == "<=") & std_values.notna() & upper.isna()
        if le_mask.any():
            upper.loc[le_mask] = std_values.loc[le_mask]

        between_mask = relations.isin({"between", "range"})
        between_ready = between_mask & std_values.notna() & std_upper.notna()
        if between_ready.any():
            between_lower = pd.Series(
                np.minimum(std_values, std_upper), index=result.index
            )
            between_upper = pd.Series(
                np.maximum(std_values, std_upper), index=result.index
            )
            lower_needed = between_ready & lower.isna()
            upper_needed = between_ready & upper.isna()
            if lower_needed.any():
                lower.loc[lower_needed] = between_lower.loc[lower_needed]
            if upper_needed.any():
                upper.loc[upper_needed] = between_upper.loc[upper_needed]

        missing_second = between_mask & (std_values.isna() | std_upper.isna())
        if missing_second.any():
            logger.warning(
                "activity_bounds_missing_range_values",
                rows=int(missing_second.sum()),
            )

    unknown_mask = relations.notna() & ~relations.isin(known_relations)
    if bounds_cfg.log_unknown_relations and unknown_mask.any():
        unique_unknown = sorted({r for r in relations[unknown_mask] if r is not None})
        logger.warning(
            "activity_bounds_unknown_relation",
            relations=unique_unknown,
            rows=int(unknown_mask.sum()),
        )

    conflict_mask = lower.notna() & upper.notna() & (lower > upper)
    if conflict_mask.any():
        logger.warning(
            "activity_bounds_conflict_swapped",
            rows=int(conflict_mask.sum()),
        )
        tmp = lower.loc[conflict_mask].copy()
        lower.loc[conflict_mask] = upper.loc[conflict_mask]
        upper.loc[conflict_mask] = tmp

    if bounds_cfg.clamp_nonnegative:
        std_types = result.get(
            "standard_type", pd.Series([None] * len(result), index=result.index)
        )
        std_units = result.get(
            "standard_units", pd.Series([None] * len(result), index=result.index)
        )
        clamp_mask = std_types.combine(std_units, _needs_nonnegative_clamp)
        clamp_series = clamp_mask.fillna(False).astype(bool)
        lower_neg = clamp_series & lower.notna() & (lower < 0)
        upper_neg = clamp_series & upper.notna() & (upper < 0)
        if lower_neg.any() or upper_neg.any():
            logger.warning(
                "activity_bounds_clamped_nonnegative",
                lower_rows=int(lower_neg.sum()),
                upper_rows=int(upper_neg.sum()),
            )
            lower.loc[lower_neg] = 0.0
            upper.loc[upper_neg] = 0.0

    digits = bounds_cfg.rounding_digits
    if digits is not None:
        lower = lower.round(digits)
        upper = upper.round(digits)

    result["lower_value"] = lower
    result["upper_value"] = upper
    return result


def _normalise_text(value: object) -> str | None:
    """Return a stripped string representation or ``None`` for blanks."""

    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    text = str(value).strip()
    return text or None


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


def _normalise_token(value: object) -> str | None:
    """Return a canonical lookup key for *value* or ``None``."""

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
    """Return a cleaned string representation of *value* or ``None``."""

    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _has_value(value: Any) -> bool:
    """Return ``True`` if *value* is considered non-empty."""

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
    """Recursively remove empty values from nested structures."""

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
    """Return concatenated text fields from *record* for keyword scans."""

    parts: list[str] = []
    for field in _COMMENT_FIELDS:
        value = record.get(field)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                parts.append(stripped)
    return " ".join(parts).lower()


def _extract_effect_features(record: Mapping[str, Any]) -> dict[str, Any]:
    """Derive allosteric effect flags from textual annotations."""

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
    """Return the first textual value from *values* that is not blank."""

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
    """Return mapped value based on *fields* or ``None`` if not found."""

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
    """Return the canonical action type for *record* using *cfg* heuristics."""

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
    """Return JSON payload and optional hash for ``activity_properties``."""

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
        triage_payload = {
            "label": _normalise_output(mapped),
            "source_field": field,
            "source_value": _normalise_output(raw_value),
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


def _normalise_mapping(mapping: Mapping[str, str]) -> dict[str, str]:
    """Return a mapping with normalised keys and stripped values."""

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
    """Return ``df`` with ``action_type`` and ``activity_properties`` columns."""

    if df.empty:
        result = df.copy()
        if action_cfg.enabled:
            result[action_cfg.column] = pd.Series(dtype="string")
        if properties_cfg.enabled:
            result[properties_cfg.column] = pd.Series(dtype="string")
        if properties_cfg.hash_column:
            result[properties_cfg.hash_column] = pd.Series(dtype="string")
        return result

    metrics_map = _normalise_mapping(action_cfg.metrics)
    functionality_map = _normalise_mapping(action_cfg.functionality)
    mechanism_map = _normalise_mapping(action_cfg.mechanism)
    triage_map = _normalise_mapping(action_cfg.triages)
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
        features = _extract_effect_features(record)
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


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    The resulting CSV places columns defined in :data:`ActivitiesSchema`
    first, preserving their declared order. Any additional fields appear
    afterwards sorted alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="activity.limit", limit=limit)
        return 1

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=expected)
        return 0

    # Configure HTTP session with the supplied User-Agent and retry policy
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv, column=cfg.activity.column, cfg=cfg.io
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "read_fail",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        # Apply the ``limit`` without materialising the entire iterator first.
        # ``itertools.islice`` allows lazy slicing; converting to ``list`` enables
        # length calculation for logging purposes.
        ids: Iterable[str] = ids_iter
        if limit is not None:
            limited = list(islice(ids_iter, limit))
            ids = limited
            logger.info("process_limit", limit=len(limited))

        enrichment_cfg = cfg.activity_enrichment
        extra_columns: list[str] = []
        action_cfg = enrichment_cfg.action_type
        if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
            extra_columns.append(action_cfg.column)
        extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}

        try:
            df = cl.get_activities(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=cfg.activity.batch_size,
                timeout=cfg.activity.timeout,
                **extra_kwargs,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "activity_fetch_failed",
                extra={"msg": str(exc)},
                error=str(exc),
                batch_size=cfg.activity.batch_size,
                timeout=cfg.activity.timeout,
            )
            return 1
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_activities(df)
        df = add_pipeline_metadata(df)
        df = compute_activity_bounds(df, cfg.activity_bounds)
        df = apply_activity_annotations(
            df,
            action_cfg=enrichment_cfg.action_type,
            properties_cfg=enrichment_cfg.activity_properties,
        )
        # Determine final column order: schema-defined columns first in their
        # declared sequence, followed by any additional columns sorted
        # alphabetically to provide deterministic output.
        schema_cols = list(ActivitiesSchema.columns)
        head = [c for c in schema_cols if c in df.columns]
        tail = sorted(c for c in df.columns if c not in schema_cols)
        col_order = head + tail
        rows_total = len(df)
        exit_code = 0
        required_cols = {
            name for name, col in ActivitiesSchema.columns.items() if col.required
        }
        optional_cols = set(ActivitiesSchema.columns) - required_cols
        missing_required = required_cols - set(df.columns)
        missing_optional = optional_cols - set(df.columns)
        if not missing_required:
            if missing_optional:
                logger.warning(
                    "optional_columns_missing",
                    columns=sorted(missing_optional),
                )
            try:
                validation_result = validate_activities(df, return_result=True)
            except SchemaErrors as exc:
                failure_path = Path(output).with_name(
                    f"{Path(output).stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path)
                logger.error(
                    "validation_failed",
                    failures=len(exc.failure_cases),
                    path=str(failure_path),
                )
                df = getattr(exc, "validated_data", df)
                exit_code = 1
            else:
                df = validation_result.data
                if not validation_result.failure_cases.empty:
                    failure_path = Path(output).with_name(
                        f"{Path(output).stem}_failure_cases.csv"
                    )
                    errors = SidecarErrors()
                    for row in validation_result.failure_cases.to_dict("records"):
                        errors.add_error(row)
                    errors.save(failure_path)
                    logger.error(
                        "validation_failed",
                        failures=len(validation_result.failure_cases),
                        path=str(failure_path),
                    )
                    exit_code = 1
        else:
            logger.warning(
                "validation_skipped",
                missing_columns=sorted(missing_required),
            )
        rows_kept = len(df)
        rows_dropped = rows_total - rows_kept
        try:
            key_cols = [c for c in ["activity_id"] if c in df.columns]
            csv_path = io.write_csv(
                df,
                output,
                cfg=cfg,
                key_cols=key_cols or None,
                col_order=col_order,
            )
            logger.info("write_done", rows=rows_kept, path=str(csv_path))
        except OSError as exc:
            logger.error(
                "write_fail",
                error=str(exc),
                path=str(output),
            )
            return 1

        stats: Stats = {
            "rows_total": rows_total,
            "rows_kept": rows_kept,
            "rows_dropped": rows_dropped,
            "output_sha256": file_sha256(csv_path),
        }
        write_meta_yaml(
            csv_path=csv_path,
            command=" ".join(sys.argv),
            config_subset=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(args.input_csv)},
            stats=stats,
            schema="ActivitiesSchema",
        )

        try:
            analyze_table_quality(df, table_name=str(output.with_suffix("")))
        except ValueError as exc:
            logger.error(
                "quality_report_failed",
                error=str(exc),
                path=str(output),
            )
            return 1
        return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read input and exit without contacting the API or writing files",
    )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    if args.limit is not None and args.limit <= 0:
        # Reject non-positive limits early to provide clear CLI feedback.
        parser.error("--limit must be a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "timeout": "activity.timeout",
                "column": "activity.column",
                "batch_size": "activity.batch_size",
                "limit": "activity.limit",
                "dry_run": "activity.dry_run",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.output_csv),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
