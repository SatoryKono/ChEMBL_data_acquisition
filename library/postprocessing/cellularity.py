"""Cellularity helpers mirroring the legacy Power Query logic."""

from __future__ import annotations

import pandas as pd

from library.pipelines.target import organism_classification
from library.schemas.targets import CELLULARITY_COLUMN_NAME

TYPE_MULTICELLULAR = "Multicellular organism"
TYPE_UNICELLULAR = "Unicellular organism"
TYPE_VIRAL = "Viral"

_VIRAL_LEGACY_LABEL = organism_classification.TYPE_VIRAL


def _normalise_label(value: object | None) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text:
        return ""
    if text.lower() in {"viruses", "virus", "viral"}:
        return TYPE_VIRAL
    if text.lower() == TYPE_UNICELLULAR.lower():
        return TYPE_UNICELLULAR
    if text.lower() == TYPE_MULTICELLULAR.lower():
        return TYPE_MULTICELLULAR
    return text


def normalise_series(series: pd.Series) -> pd.Series:
    """Return ``series`` with canonical cellularity labels."""

    return series.map(_normalise_label).astype("string")


def ensure_cellularity(
    frame: pd.DataFrame,
    *,
    genus_col: str = organism_classification.DEFAULT_GENUS_COLUMN,
    superkingdom_col: str = organism_classification.DEFAULT_SUPERKINGDOM_COLUMN,
    phylum_col: str = organism_classification.DEFAULT_PHYLUM_COLUMN,
    class_col: str = organism_classification.DEFAULT_CLASS_COLUMN,
    output_col: str = CELLULARITY_COLUMN_NAME,
) -> pd.DataFrame:
    """Return ``frame`` with the canonical cellularity column present."""

    result = frame.copy()
    if output_col not in result.columns:
        result = organism_classification.add_cellularity_smart(
            result,
            genus_col=genus_col,
            superkingdom_col=superkingdom_col,
            phylum_col=phylum_col,
            class_col=class_col,
            output_col=output_col,
        )
    result[output_col] = normalise_series(result[output_col])
    return result


def unicellular_flag(series: pd.Series) -> pd.Series:
    """Return the boolean unicellular flag derived from ``series`` labels."""

    lowered = series.astype("string").str.strip().str.lower()
    mask = lowered.isin(
        {TYPE_UNICELLULAR.lower(), TYPE_VIRAL.lower(), _VIRAL_LEGACY_LABEL.lower()}
    )
    return mask.astype("boolean")
