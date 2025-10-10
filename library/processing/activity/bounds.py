"""Utilities for computing activity measurement bounds."""

from __future__ import annotations

import re
from collections.abc import Iterable

import numpy as np
import pandas as pd

from ...common.log import logger
from ...common.pandas_utils import merge_series_prefer_left
from ...config import ActivityBoundsCfg

__all__ = ["compute_activity_bounds"]

_COVERAGE_WARN_THRESHOLD = 0.95

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


def _swap_conflicts(lower: pd.Series, upper: pd.Series) -> None:
    conflict_mask = lower.notna() & upper.notna() & (lower > upper)
    if conflict_mask.any():
        logger.warning(
            "activity_bounds_conflict_swapped",
            rows=int(conflict_mask.sum()),
        )
        tmp = lower.loc[conflict_mask].copy()
        lower.loc[conflict_mask] = upper.loc[conflict_mask]
        upper.loc[conflict_mask] = tmp


def _clamp_nonnegative(
    df: pd.DataFrame, lower: pd.Series, upper: pd.Series
) -> tuple[pd.Series, pd.Series]:
    std_types = df.get("standard_type", pd.Series([None] * len(df), index=df.index))
    std_units = df.get("standard_units", pd.Series([None] * len(df), index=df.index))
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
    return lower, upper


def _log_bounds_coverage(
    std_lower: pd.Series,
    std_upper: pd.Series,
    lower: pd.Series,
    upper: pd.Series,
) -> None:
    total = len(lower)
    if total == 0:
        return

    def _coverage(series: pd.Series) -> tuple[int, float]:
        count = int(series.notna().sum())
        ratio = count / total if total else 0.0
        return count, ratio

    std_lower_count, std_lower_ratio = _coverage(std_lower)
    std_upper_count, std_upper_ratio = _coverage(std_upper)
    lower_count, lower_ratio = _coverage(lower)
    upper_count, upper_ratio = _coverage(upper)

    logger.info(
        "activity_bounds_coverage",
        rows=total,
        standard_lower_count=std_lower_count,
        standard_lower_pct=round(std_lower_ratio * 100, 2),
        standard_upper_count=std_upper_count,
        standard_upper_pct=round(std_upper_ratio * 100, 2),
        lower_count=lower_count,
        lower_pct=round(lower_ratio * 100, 2),
        upper_count=upper_count,
        upper_pct=round(upper_ratio * 100, 2),
    )

    for column, input_ratio, output_ratio in (
        ("lower_value", std_lower_ratio, lower_ratio),
        ("upper_value", std_upper_ratio, upper_ratio),
    ):
        if (
            input_ratio >= _COVERAGE_WARN_THRESHOLD
            and output_ratio < _COVERAGE_WARN_THRESHOLD
        ):
            logger.warning(
                "activity_bounds_low_output_coverage",
                column=column,
                rows=total,
                input_pct=round(input_ratio * 100, 2),
                output_pct=round(output_ratio * 100, 2),
            )


def _normalize_relations(df: pd.DataFrame) -> pd.Series:
    return _relation_series(df).map(_normalize_relation)


def _log_unknown_relations(relations: pd.Series, known: Iterable[str]) -> None:
    known_set = set(known)
    unknown_mask = relations.notna() & ~relations.isin(known_set)
    if unknown_mask.any():
        unique_unknown = sorted({r for r in relations[unknown_mask] if r is not None})
        logger.warning(
            "activity_bounds_unknown_relation",
            relations=unique_unknown,
            rows=int(unknown_mask.sum()),
        )


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

    lower = merge_series_prefer_left(lower, std_lower)
    upper = merge_series_prefer_left(upper, std_upper)

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

    relations = _normalize_relations(result)
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

    if bounds_cfg.log_unknown_relations:
        _log_unknown_relations(relations, known_relations)

    _swap_conflicts(lower, upper)

    if bounds_cfg.clamp_nonnegative:
        lower, upper = _clamp_nonnegative(result, lower, upper)

    digits = bounds_cfg.rounding_digits
    if digits is not None:
        lower = lower.round(digits)
        upper = upper.round(digits)

    result["lower_value"] = lower
    result["upper_value"] = upper
    if "activity_id" in result.columns and pd.api.types.is_string_dtype(
        result["activity_id"].dtype
    ):
        result["activity_id"] = result["activity_id"].astype("object")
    _log_bounds_coverage(
        std_lower, std_upper, result["lower_value"], result["upper_value"]
    )
    return result
