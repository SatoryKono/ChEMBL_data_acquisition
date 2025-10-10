"""Utility helpers shared across target post-processing steps."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any
from xml.etree.ElementTree import Element

import pandas as pd


def normalize_text(value: Any) -> str:
    """Replicate the Power Query ``NormalizeText`` helper."""

    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    if value is pd.NA:
        return ""
    text = str(value)
    return text.strip().lower()


def first_element_text(nodes: Iterable[Element], element_name: str) -> str | None:
    """Return the text of the first XML node whose tag matches ``element_name``."""

    for node in nodes:
        if node.tag != element_name:
            continue
        text = node.text
        if text is None:
            return None
        return text
    return None


__all__ = ["normalize_text", "first_element_text"]
