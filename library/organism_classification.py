"""Utilities for classifying organism cellularity.

This module provides helpers to derive a ``unicellular_organism`` flag
from basic taxonomic annotations available in ChEMBL dictionaries.
"""

from __future__ import annotations

from typing import Iterable

import pandas as pd

# Known taxonomic groups that are predominantly unicellular.  The values are
# stored in lower case to simplify string normalisation via ``str.casefold``.
_UNICELLULAR_SUPERKINGDOMS: set[str] = {
    "bacteria",
    "archaea",
    "viruses",
}

_UNICELLULAR_EUKARYOTIC_PHYLA: set[str] = {
    "sar",
    "discoba",
}

_UNICELLULAR_EUKARYOTIC_GENERA: set[str] = {
    "candida",
    "chlamydomonas",
    "crithidia",
    "cryptosporidium",
    "eimeria",
    "leishmania",
    "malassezia",
    "plasmodium (laverania)",
    "plasmodium (plasmodium)",
    "pneumocystis",
    "schizotrypanum",
    "toxoplasma",
    "trypanosoma",
}

_VIRUS_TAXON_IDS: set[str] = {
    "10246",
    "11082",
    "11676",
    "11709",
    "11908",
    "11926",
    "169066",
    "211044",
    "266827",
    "2697049",
    "3052230",
    "31647",
    "343983",
    "365124",
    "382835",
    "64320",
    "648213",
    "694009",
}


def _normalise(series: pd.Series) -> pd.Series:
    """Strip and case-fold string values for robust comparisons."""

    return (
        series.astype("string")
        .str.strip()
        .str.replace("\\s+", " ", regex=True)
        .str.casefold()
    )


def _any_in(values: pd.Series, options: Iterable[str]) -> pd.Series:
    """Check whether entries belong to a predefined set."""

    return values.isin(set(options))


def add_cellularity_smart(
    df: pd.DataFrame,
    *,
    genus_col: str,
    superkingdom_col: str,
    phylum_col: str,
    taxon_id_col: str | None = None,
    output_col: str = "unicellular_organism",
) -> pd.DataFrame:
    """Annotate unicellular organisms based on taxonomic context."""

    genus = _normalise(df.get(genus_col, pd.Series(dtype="string")))
    superkingdom = _normalise(df.get(superkingdom_col, pd.Series(dtype="string")))
    phylum = _normalise(df.get(phylum_col, pd.Series(dtype="string")))
    if taxon_id_col:
        taxon_id = _normalise(df.get(taxon_id_col, pd.Series(dtype="string")))
    else:
        taxon_id = pd.Series(pd.NA, index=df.index, dtype="string")

    unicellular = pd.Series(False, index=df.index, dtype="boolean")

    unicellular |= _any_in(superkingdom, _UNICELLULAR_SUPERKINGDOMS)
    unicellular |= _any_in(taxon_id, _VIRUS_TAXON_IDS)

    eukaryotic = superkingdom == "eukaryota"
    unicellular |= eukaryotic & _any_in(phylum, _UNICELLULAR_EUKARYOTIC_PHYLA)
    unicellular |= eukaryotic & _any_in(genus, _UNICELLULAR_EUKARYOTIC_GENERA)

    df[output_col] = unicellular.fillna(False)
    return df
