"""Helper utilities mirroring the Power Query helpers used in the organism
post-processing workbook."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any
from xml.etree import ElementTree as ET

__all__ = ["normalize_text", "first_element_text", "distinct_preserve_order"]


def normalize_text(value: Any) -> str:
    """Return lower-cased text stripped of surrounding whitespace.

    This is a direct translation of the Power Query ``NormalizeText`` helper.
    ``null`` inputs in M map to ``None`` in Python and become empty strings.
    ``Text.Lower`` and ``Text.Trim`` are applied sequentially otherwise.
    """

    if value is None:
        return ""
    text = str(value)
    return text.strip().lower()


def first_element_text(nodes: Iterable[ET.Element], element_name: str) -> str | None:
    """Return the text of the first XML element matching ``element_name``."""

    for element in nodes:
        if element.tag == element_name:
            text = element.text
            if text is not None:
                return text
            return None
    return None


def distinct_preserve_order(items: Iterable[Any]) -> list[Any]:
    """Return values from ``items`` preserving the first occurrence order."""

    seen: list[Any] = []
    for item in items:
        if item not in seen:
            seen.append(item)
    return seen
