"""Heuristics for inferring organism cellularity from taxonomy data."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

__all__ = ["add_cellularity_smart"]

TYPE_MULTICELLULAR = "Multicellular organism"
TYPE_UNICELLULAR = "Unicellular organism"
TYPE_VIRAL = "Viruses"

_VIRAL_SUPERKINGDOMS: frozenset[str] = frozenset({"viruses"})
_UNICELLULAR_SUPERKINGDOMS: frozenset[str] = frozenset({"bacteria", "archaea"})
_UNICELLULAR_PHYLUM: frozenset[str] = frozenset(
    {
        "apicomplexa",
        "amoebozoa",
        "ciliophora",
        "chlorophyta",
        "euglenozoa",
        "myzozoa",
        "metamonada",
        "microsporidia",
    }
)
_UNICELLULAR_CLASSES: frozenset[str] = frozenset(
    {
        "aconoidasida",
        "conoidasida",
        "kinetoplastea",
        "saccharomycetes",
        "pneumocystidomycetes",
        "malasseziomycetes",
        "chlorophyceae",
    }
)
_UNICELLULAR_GENUS: frozenset[str] = frozenset(
    {
        "candida",
        "malassezia",
        "pneumocystis",
        "chlamydomonas",
    }
)


def _normalise(value: object | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    lowered = text.lower()
    if lowered in {"nan", "none", "-", "na"}:
        return ""
    return lowered


def _tokens(value: str) -> Iterable[str]:
    if not value:
        return ()
    return [token.strip() for token in value.replace(";", "|").split("|") if token.strip()]


def _classify_row(
    genus: str,
    superkingdom: str,
    phylum: str,
    klass: str,
) -> str:
    if superkingdom in _VIRAL_SUPERKINGDOMS:
        return TYPE_VIRAL
    if superkingdom in _UNICELLULAR_SUPERKINGDOMS:
        return TYPE_UNICELLULAR

    if any(token in _UNICELLULAR_PHYLUM for token in _tokens(phylum)):
        return TYPE_UNICELLULAR
    if any(token in _UNICELLULAR_CLASSES for token in _tokens(klass)):
        return TYPE_UNICELLULAR
    if genus in _UNICELLULAR_GENUS:
        return TYPE_UNICELLULAR
    return TYPE_MULTICELLULAR


def add_cellularity_smart(
    df: pd.DataFrame,
    *,
    genus_col: str = "genus",
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str = "lineage_class",
    output_col: str = "type",
) -> pd.DataFrame:
    """Return ``df`` with an inferred organism cellularity column.

    The classification relies on UniProt taxonomy fields. Viruses are detected
    from the superkingdom, bacteria and archaea default to unicellular and a set
    of unicellular eukaryote phyla/classes cover protozoa, yeasts and algae.

    Parameters
    ----------
    df:
        Input table containing taxonomy information.
    genus_col, superkingdom_col, phylum_col, class_col:
        Column names holding genus and taxonomy lineage values.
    output_col:
        Name of the created column. Defaults to ``"type"``.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with the ``output_col`` column appended. The column uses
        string dtype.
    """

    result = df.copy()
    if result.empty:
        result[output_col] = pd.Series(dtype="string")
        return result

    def _infer(row: pd.Series) -> str:
        genus = _normalise(row.get(genus_col))
        superkingdom = _normalise(row.get(superkingdom_col))
        phylum = _normalise(row.get(phylum_col))
        klass = _normalise(row.get(class_col))
        return _classify_row(genus, superkingdom, phylum, klass)

    result[output_col] = result.apply(_infer, axis=1).astype("string")
    return result
