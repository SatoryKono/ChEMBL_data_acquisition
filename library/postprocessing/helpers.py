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
from pathlib import Path
from typing import Iterable, Sequence
from xml.etree import ElementTree

import pandas as pd

from library.common.csv_utils import write_csv_deterministic

# Accepted encodings – these match the Power Query ``Binary.Decompress`` fallbacks
# that were historically used when loading the aggregated targets table.
ENCODING_FALLBACKS: tuple[str, ...] = ("utf-8", "utf-8-sig", "cp1252")


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
    try:
        if pd.isna(value):  # type: ignore[arg-type]
            return ""
    except TypeError:
        pass
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
            return pd.read_csv(path, sep=sep, encoding=encoding, dtype="string")
        except Exception as exc:  # pragma: no cover - defensive guard
            errors.append(exc)
    if not errors:  # pragma: no cover - defensive guard
        raise RuntimeError(f"Unable to read CSV at {path!s}")
    # Re-raise the last error so the traceback points to the failing codec.
    raise errors[-1]


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


def fill_missing(frame: pd.DataFrame, columns: Iterable[str], fill_value: str = "-") -> pd.DataFrame:
    """Fill ``columns`` missing from ``frame`` with ``fill_value``."""

    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = fill_value
    return result

