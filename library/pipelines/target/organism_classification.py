"""Utilities for classifying organism cellularity.

This module provides helpers to derive a ``unicellular_organism`` flag
from basic taxonomic annotations available in ChEMBL dictionaries.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, cast

import pandas as pd

__all__ = [
    "OrganismClassificationRules",
    "normalize",
    "classify_by_lineage",
    "classify_record",
    "add_cellularity",
    "add_cellularity_smart",
    "TYPE_MULTICELLULAR",
    "TYPE_UNICELLULAR",
    "TYPE_VIRAL",
]

# ===== Text labels =====
TYPE_MULTICELLULAR = "Multicellular organism"
TYPE_UNICELLULAR = "Unicellular organism"
TYPE_VIRAL = "Viruses"

# ===== Defaults =====
DEFAULT_GENUS_COLUMN = "genus"
DEFAULT_SUPERKINGDOM_COLUMN = "lineage_superkingdom"
DEFAULT_PHYLUM_COLUMN = "lineage_phylum"
DEFAULT_CLASS_COLUMN = "lineage_class"
DEFAULT_OUTPUT_COLUMN = "type"
DEFAULT_TAXON_ID_COLUMN = "taxon_id"

DEFAULT_NULL_LITERALS: frozenset[str] = frozenset({"nan", "none", "-", "na"})
TAXONOMY_SEPARATORS: Sequence[str] = ("|", ";")
EMPTY_VALUE = ""

# ===== Rule lists =====
DEFAULT_VIRAL_SUPERKINGDOMS: frozenset[str] = frozenset({"viruses"})
DEFAULT_UNICELLULAR_SUPERKINGDOMS: frozenset[str] = frozenset({"bacteria", "archaea"})
DEFAULT_UNICELLULAR_PHYLA: frozenset[str] = frozenset(
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
DEFAULT_UNICELLULAR_CLASSES: frozenset[str] = frozenset(
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
DEFAULT_UNICELLULAR_GENERA: frozenset[str] = frozenset(
    {
        "candida",
        "malassezia",
        "pneumocystis",
        "chlamydomonas",
    }
)


@dataclass(frozen=True)
class OrganismClassificationRules:
    """Configuration container for cellularity heuristics."""

    viral_superkingdoms: frozenset[str] = field(
        default_factory=lambda: DEFAULT_VIRAL_SUPERKINGDOMS
    )
    unicellular_superkingdoms: frozenset[str] = field(
        default_factory=lambda: DEFAULT_UNICELLULAR_SUPERKINGDOMS
    )
    unicellular_phyla: frozenset[str] = field(
        default_factory=lambda: DEFAULT_UNICELLULAR_PHYLA
    )
    unicellular_classes: frozenset[str] = field(
        default_factory=lambda: DEFAULT_UNICELLULAR_CLASSES
    )
    unicellular_genera: frozenset[str] = field(
        default_factory=lambda: DEFAULT_UNICELLULAR_GENERA
    )
    labels: Sequence[str] = field(
        default_factory=lambda: (TYPE_MULTICELLULAR, TYPE_UNICELLULAR, TYPE_VIRAL)
    )
    null_literals: frozenset[str] = field(default_factory=lambda: DEFAULT_NULL_LITERALS)

    def label_multicellular(self) -> str:
        return self.labels[0]

    def label_unicellular(self) -> str:
        return self.labels[1]

    def label_viral(self) -> str:
        return self.labels[2]


DEFAULT_RULES = OrganismClassificationRules()


def normalize(
    value: object | None, *, rules: OrganismClassificationRules = DEFAULT_RULES
) -> str:
    """Return a lowercase normalised taxonomy token."""

    if isinstance(value, pd.Series):
        raise TypeError("normalize expects a scalar value")
    if value is None:
        return EMPTY_VALUE
    try:
        if pd.isna(cast(Any, value)):
            return EMPTY_VALUE
    except TypeError:
        pass

    text = str(value).strip()
    if not text:
        return EMPTY_VALUE
    lowered = text.lower()
    if lowered in rules.null_literals:
        return EMPTY_VALUE
    return lowered


def _split_taxonomy(value: str) -> Iterable[str]:
    if not value:
        return ()

    tokens = [value]
    for separator in TAXONOMY_SEPARATORS:
        tokens = [token for chunk in tokens for token in chunk.split(separator)]
    return (token.strip() for token in tokens if token.strip())


def classify_by_lineage(
    *,
    genus: str,
    superkingdom: str,
    phylum: str,
    klass: str,
    rules: OrganismClassificationRules = DEFAULT_RULES,
) -> str:
    """Classify an organism using lineage information."""

    if superkingdom in rules.viral_superkingdoms:
        return rules.label_viral()
    if superkingdom in rules.unicellular_superkingdoms:
        return rules.label_unicellular()

    if any(token in rules.unicellular_phyla for token in _split_taxonomy(phylum)):
        return rules.label_unicellular()
    if any(token in rules.unicellular_classes for token in _split_taxonomy(klass)):
        return rules.label_unicellular()
    if genus in rules.unicellular_genera:
        return rules.label_unicellular()

    return rules.label_multicellular()


def classify_record(
    record: Mapping[str, object],
    *,
    genus_column: str = DEFAULT_GENUS_COLUMN,
    superkingdom_column: str = DEFAULT_SUPERKINGDOM_COLUMN,
    phylum_column: str = DEFAULT_PHYLUM_COLUMN,
    class_column: str = DEFAULT_CLASS_COLUMN,
    rules: OrganismClassificationRules = DEFAULT_RULES,
) -> str:
    """Infer the cellularity label for a mapping of lineage fields."""

    genus = normalize(record.get(genus_column), rules=rules)
    superkingdom = normalize(record.get(superkingdom_column), rules=rules)
    phylum = normalize(record.get(phylum_column), rules=rules)
    klass = normalize(record.get(class_column), rules=rules)
    return classify_by_lineage(
        genus=genus,
        superkingdom=superkingdom,
        phylum=phylum,
        klass=klass,
        rules=rules,
    )


def add_cellularity(
    df: pd.DataFrame,
    *,
    genus_column: str = DEFAULT_GENUS_COLUMN,
    superkingdom_column: str = DEFAULT_SUPERKINGDOM_COLUMN,
    phylum_column: str = DEFAULT_PHYLUM_COLUMN,
    class_column: str = DEFAULT_CLASS_COLUMN,
    output_column: str = DEFAULT_OUTPUT_COLUMN,
    rules: OrganismClassificationRules = DEFAULT_RULES,
) -> pd.DataFrame:
    """Return ``df`` with the inferred cellularity column appended."""

    result = df.copy()
    if result.empty:
        result[output_column] = pd.Series(dtype="string")
        return result

    def _infer(row: pd.Series) -> str:
        return classify_record(
            row,
            genus_column=genus_column,
            superkingdom_column=superkingdom_column,
            phylum_column=phylum_column,
            class_column=class_column,
            rules=rules,
        )

    result[output_column] = result.apply(_infer, axis=1).astype("string")
    return result


def add_cellularity_smart(
    df: pd.DataFrame,
    *,
    genus_col: str = DEFAULT_GENUS_COLUMN,
    superkingdom_col: str = DEFAULT_SUPERKINGDOM_COLUMN,
    phylum_col: str = DEFAULT_PHYLUM_COLUMN,
    class_col: str = DEFAULT_CLASS_COLUMN,
    output_col: str = DEFAULT_OUTPUT_COLUMN,
    taxon_id_col: str = DEFAULT_TAXON_ID_COLUMN,
    rules: OrganismClassificationRules = DEFAULT_RULES,
) -> pd.DataFrame:
    """Backward-compatible wrapper around :func:`add_cellularity`.

    Parameters
    ----------
    df:
        Input frame expected to contain the lineage columns used for
        classification.
    genus_col, superkingdom_col, phylum_col, class_col:
        Column names storing taxonomy annotations.
    output_col:
        Column receiving the inferred cellularity labels.
    taxon_id_col:
        Retained for backwards compatibility with legacy call sites. The
        identifier is not required for lineage-based inference and therefore
        ignored.
    rules:
        Rule set controlling the inference heuristics.
    """

    _ = taxon_id_col  # Intentionally ignored; maintained for compatibility.

    return add_cellularity(
        df,
        genus_column=genus_col,
        superkingdom_column=superkingdom_col,
        phylum_column=phylum_col,
        class_column=class_col,
        output_column=output_col,
        rules=rules,
    )
