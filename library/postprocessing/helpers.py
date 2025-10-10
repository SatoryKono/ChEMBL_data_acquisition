"""Common helpers for target post-processing.

The original pipeline relied on a Power Query workbook where a number of text
cleaning and XML parsing helpers were reused across multiple steps.  The
implementation below mirrors those routines so that the Python port remains
byte-for-byte compatible with the canonical Single Source of Truth (SSoT)
exports.  The helpers are intentionally small and focused: they deal with text
normalisation, best-effort XML parsing and deterministic CSV emission.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TypeGuard
from xml.etree import ElementTree

import numpy as np
import pandas as pd
from pandas._typing import Scalar
from pandas.errors import ParserError

from library.common.csv_utils import write_csv_deterministic

# Accepted encodings – these match the Power Query ``Binary.Decompress`` fallbacks
# that were historically used when loading the aggregated targets table.
# ``latin-1`` (and its canonical alias ``iso-8859-1``) are included because a
# handful of legacy dictionary exports produced by Power Query use extended
# control characters (e.g. ``0x81``) that are not representable in ``cp1252``.
# The ISO-8859-1 codec accepts these bytes and preserves them verbatim so the
# downstream normalisation logic can continue operating deterministically even
# on systems where ``latin-1`` is not registered as a codec (most notably some
# Windows Python distributions).
ENCODING_FALLBACKS: tuple[str, ...] = (
    "utf-8",
    "utf-8-sig",
    "cp1252",
    "latin-1",
    "iso-8859-1",
)
CSV_SEPARATORS: tuple[str, ...] = (",", "\t", ";")


def _is_scalar_like(value: object) -> TypeGuard[Scalar]:
    """Return ``True`` when ``value`` is a pandas scalar."""

    return bool(pd.api.types.is_scalar(value))


def is_missing_scalar(value: object | None) -> bool:
    """Return ``True`` when ``value`` should be treated as a missing scalar."""

    if value is None:
        return True
    if _is_scalar_like(value):
        return bool(pd.isna(value))
    return False


def normalise_export_basename(path: Path) -> str:
    """Return a deterministic basename for post-processed exports.

    The helper mirrors the naming conventions used by the legacy Power Query
    workbook: temporary suffixes (``.tmp``) are stripped, leading dots are
    removed so that subsequent ``Path.with_name`` calls do not introduce
    accidental hidden files, and the remaining stem is returned unchanged.
    """

    name = path.name
    while name.startswith("."):
        name = name[1:]
    tmp_suffix = ".tmp"
    while name.lower().endswith(tmp_suffix):
        name = name[: -len(tmp_suffix)]
    stem, ext = os.path.splitext(name)
    ext_lower = ext.lower()
    for suffix in ("_normalized", "_normalised"):
        if ext_lower.endswith(suffix):
            ext = ext[: -len(suffix)] or ""
            if ext == ".":
                ext = ""
            ext_lower = ext.lower()
            break
    lowered_stem = stem.lower()
    for suffix in ("_normalized", "_normalised"):
        if lowered_stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    name = f"{stem}{ext or ''}"
    if not name:
        raise ValueError(f"Unable to derive export basename from {path!s}")
    return name


def normalise_text(value: object | None) -> str:
    """Return a trimmed textual representation of ``value``.

    ``None`` values, NaNs and empty strings are mapped to ``""`` mirroring the
    behaviour of the Power Query ``ToText`` helper.
    """

    if value is None:
        return ""
    if isinstance(value, str):
        text = value.strip()
        return text
    if is_missing_scalar(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"none", "nan"} else text


def ensure_string_columns(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    """Return ``frame`` with the specified ``columns`` coerced to ``string`` dtype."""

    result = frame.copy()
    for column in columns:
        if column in result.columns:
            result[column] = result[column].astype("string")
    return result


def read_csv_with_fallbacks(
    path: Path | str,
    *,
    sep: str = ",",
    encodings: Sequence[str] = ENCODING_FALLBACKS,
) -> pd.DataFrame:
    """Load a CSV file trying the known Power Query encoding fallbacks."""

    errors: list[Exception] = []
    for encoding in encodings:
        try:
            frame = pd.read_csv(path, sep=sep, encoding=encoding, dtype="string")
            # Normalise column headers so UTF-8 with BOM exports behave consistently.
            frame = frame.rename(columns=lambda c: str(c).lstrip("\ufeff"))
            return frame
        except Exception as exc:  # pragma: no cover - defensive guard
            errors.append(exc)
    if not errors:  # pragma: no cover - defensive guard
        raise RuntimeError(f"Unable to read CSV at {path!s}")
    # Re-raise the last error so the traceback points to the failing codec.
    raise errors[-1]


def _normalise_encoding_candidates(
    encoding: str | Sequence[str] | None,
) -> list[str]:
    candidates: list[str] = []
    seen: set[str] = set()

    def _append(value: str) -> None:
        key = value.lower()
        if key not in seen:
            seen.add(key)
            candidates.append(value)

    if encoding is None:
        pass
    elif isinstance(encoding, str):
        if encoding:
            _append(encoding)
    else:
        for item in encoding:
            if item:
                _append(item)

    for fallback in ENCODING_FALLBACKS:
        _append(fallback)

    if not candidates:
        candidates.append("utf-8")
    return candidates


def _normalise_separator_candidates(separators: Sequence[str] | None) -> list[str]:
    if separators is None:
        return list(CSV_SEPARATORS)
    seen: set[str] = set()
    resolved: list[str] = []
    for separator in (*separators, *CSV_SEPARATORS):
        if not separator:
            continue
        if separator in seen:
            continue
        seen.add(separator)
        resolved.append(separator)
    return resolved or [","]


def _dtype_for_read(
    schema: Mapping[str, str] | None,
) -> dict[str, str | pd.api.extensions.ExtensionDtype]:
    if not schema:
        return {}
    mapping: dict[str, str | pd.api.extensions.ExtensionDtype] = {}
    for column, declared in schema.items():
        kind = str(declared).strip().lower()
        if kind in {"text", "logical"}:
            mapping[column] = "string"
        elif kind == "int64":
            mapping[column] = pd.Int64Dtype()
        elif kind == "float64":
            mapping[column] = pd.Float64Dtype()
        else:
            mapping[column] = declared
    return mapping


def _coerce_logical_series(series: pd.Series) -> pd.Series:
    if series.dtype == "boolean":
        return series

    def convert(value: object) -> object:
        if is_missing_scalar(value):
            return pd.NA
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if isinstance(value, (int, np.integer)):
            if value in (0, 1):
                return bool(value)
        if isinstance(value, float):
            if value in (0.0, 1.0):
                return bool(int(value))
        text = str(value).strip()
        if not text or text == "-":
            return pd.NA
        lowered = text.lower()
        if lowered in {"true", "t", "1", "y", "yes"}:
            return True
        if lowered in {"false", "f", "0", "n", "no"}:
            return False
        raise ValueError(f"Cannot coerce value {value!r} to Logical")

    mapped = series.map(convert)
    return mapped.astype("boolean")


def _coerce_int64_series(series: pd.Series) -> pd.Series:
    if series.dtype == "Int64":
        return series
    if series.empty:
        return series.astype("Int64")
    if series.dtype == "float64":
        numeric = series.astype("float64")
    else:
        as_text = series.astype("string").str.strip()
        numeric = as_text.replace({"": pd.NA, "-": pd.NA})
    result = pd.to_numeric(numeric, errors="raise")
    return result.astype("Int64")


def _coerce_float64_series(series: pd.Series) -> pd.Series:
    if series.dtype == "Float64":
        return series
    if series.empty:
        return series.astype("Float64")
    numeric = pd.to_numeric(series, errors="coerce")
    return pd.Series(numeric, index=series.index, dtype="Float64")


def coerce_types(df: pd.DataFrame, schema: Mapping[str, str] | None) -> pd.DataFrame:
    """Return ``df`` with columns coerced according to ``schema``."""

    if not schema:
        return df.copy()

    coerced = df.copy()
    for column, declared in schema.items():
        if column not in coerced.columns:
            continue
        kind = str(declared).strip().lower()
        if kind == "text":
            coerced[column] = coerced[column].astype("string")
        elif kind == "logical":
            coerced[column] = _coerce_logical_series(coerced[column])
        elif kind == "int64":
            coerced[column] = _coerce_int64_series(coerced[column])
        elif kind == "float64":
            coerced[column] = _coerce_float64_series(coerced[column])
        else:
            raise ValueError(
                f"Unsupported schema type {declared!r} for column '{column}'"
            )
    return coerced


def read_csv_strict(
    path: Path | str,
    *,
    encoding: str | Sequence[str] | None,
    dtype_map: Mapping[str, str] | None,
    na_values: Sequence[str] | None,
    separators: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Read ``path`` enforcing workbook delimiter, NA and encoding conventions."""

    encodings = _normalise_encoding_candidates(encoding)
    resolved_separators = _normalise_separator_candidates(separators)
    dtype = _dtype_for_read(dtype_map)
    expected_columns = set(dtype_map.keys()) if dtype_map else set()
    na_tokens = list(na_values or []) or None

    last_error: Exception | None = None
    for candidate_encoding in encodings:
        for separator in resolved_separators:
            try:
                frame = pd.read_csv(
                    path,
                    sep=separator,
                    encoding=candidate_encoding,
                    dtype=dtype or None,
                    keep_default_na=False,
                    na_values=na_tokens,
                )
            except UnicodeDecodeError as exc:
                last_error = exc
                break
            except (ParserError, ValueError) as exc:
                last_error = exc
                continue
            if expected_columns and expected_columns.isdisjoint(frame.columns):
                last_error = ValueError(
                    "Detected separator did not yield required columns; trying next candidate."
                )
                continue
            frame = frame.rename(columns=lambda c: str(c).lstrip("\ufeff"))
            if dtype_map:
                frame = coerce_types(frame, dtype_map)
            return frame
    if last_error is not None:
        raise last_error
    raise RuntimeError(f"Unable to read CSV at {path!s}")


def stable_sort(df: pd.DataFrame, keys: Sequence[str]) -> pd.DataFrame:
    """Return ``df`` sorted by ``keys`` using a stable mergesort."""

    if not keys:
        return df.copy()
    return df.sort_values(by=list(keys), kind="mergesort")


def extract_xml_texts(value: str, xpath: str) -> list[str]:
    """Return the text nodes matching ``xpath`` from an XML fragment.

    The helper is intentionally permissive: parsing failures yield an empty list
    to mirror the defensive ``try ... otherwise { { "Value" = "" } }`` guards in
    the original M workbook.
    """

    text = normalise_text(value)
    if not text:
        return []
    try:
        root = ElementTree.fromstring(text)
    except ElementTree.ParseError:
        return []
    return [normalise_text(node.text) for node in root.findall(xpath)]


def write_csv(frame: pd.DataFrame, path: Path, *, columns: Sequence[str]) -> Path:
    """Serialise ``frame`` to ``path`` using deterministic ordering."""

    return write_csv_deterministic(
        frame.copy(),
        path,
        col_order=columns,
        key_cols=[columns[0]] if columns else [],
        encoding="utf-8",
        sep=",",
        cfg=None,
    )


def sort_power_query(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Return ``frame`` sorted by ``columns`` using Power Query compatible rules."""

    if frame.empty:
        return frame.copy()

    seen: list[str] = []
    for column in columns:
        if column in frame.columns and column not in seen:
            seen.append(column)

    if not seen:
        return frame.copy()

    return frame.sort_values(by=seen, kind="mergesort", na_position="last").reset_index(
        drop=True
    )


def fill_missing(
    frame: pd.DataFrame, columns: Iterable[str], fill_value: str = "-"
) -> pd.DataFrame:
    """Fill ``columns`` missing from ``frame`` with ``fill_value``."""

    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = fill_value
    return result
