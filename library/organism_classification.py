"""Classification helpers for organism cellularity labels."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import pandas as pd

# ---------------------------------------------------------------------------
# Text labels exposed by the original Power Query implementation.  They are
# defined as module-level constants to avoid magic strings in the public API
# and tests.
CELLULARITY_MULTICELLULAR = "Multicellular organism"
CELLULARITY_UNICELLULAR = "Unicellular organism"
CELLULARITY_VIRUS = "Viruses"


# ---------------------------------------------------------------------------
# Taxonomy buckets used by the legacy rules.  The values reflect the curated
# mapping delivered with the project (``dictionary/_Target/organism.csv``) and
# are intentionally conservative so that class-based overrides only trigger for
# taxa known to require them.
VIRUS_SUPERKINGDOMS = frozenset({"Viruses"})
VIRUS_PHYLA = frozenset({"Riboviria", "Varidnaviria"})

UNICELLULAR_SUPERKINGDOMS = frozenset({"Bacteria", "Archaea"})
UNICELLULAR_PHYLA = frozenset({"Discoba", "Sar"})
UNICELLULAR_CLASSES = frozenset(
    {
        "Chlorophyceae",  # unicellular green algae (e.g. Chlamydomonas)
        "Saccharomycetes",  # yeasts such as Candida
        "Pneumocystidomycetes",  # Pneumocystis jirovecii
        "Malasseziomycetes",  # lipid-dependent yeasts
    }
)

MULTICELLULAR_PHYLA = frozenset(
    {"Chordata", "Ecdysozoa", "Spiralia", "Echinodermata", "Viridiplantae"}
)
MULTICELLULAR_CLASSES = frozenset({"Agaricomycetes"})  # filamentous fungi


@dataclass(frozen=True)
class OrganismClassificationRules:
    """Container describing how taxonomy values map to cellularity labels."""

    label_multicellular: str = CELLULARITY_MULTICELLULAR
    label_unicellular: str = CELLULARITY_UNICELLULAR
    label_viruses: str = CELLULARITY_VIRUS
    virus_superkingdoms: frozenset[str] = VIRUS_SUPERKINGDOMS
    virus_phyla: frozenset[str] = VIRUS_PHYLA
    unicellular_superkingdoms: frozenset[str] = UNICELLULAR_SUPERKINGDOMS
    unicellular_phyla: frozenset[str] = UNICELLULAR_PHYLA
    unicellular_classes: frozenset[str] = UNICELLULAR_CLASSES
    multicellular_phyla: frozenset[str] = MULTICELLULAR_PHYLA
    multicellular_classes: frozenset[str] = MULTICELLULAR_CLASSES


DEFAULT_RULES = OrganismClassificationRules()


def normalize(value: object | None) -> str | None:
    """Convert raw taxonomy values to canonical text.

    ``None``, ``NaN`` and empty strings are normalised to ``None`` so callers
    can distinguish them from meaningful taxonomy entries.
    """

    if value is None:
        return None
    if isinstance(value, str):
        text = value.strip()
        if not text:
            return None
        if text.lower() in {"nan", "none", "null", "<na>"}:
            return None
        return text
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.lower() in {"nan", "none", "null", "<na>"}:
        return None
    return text


def classify_by_lineage(
    superkingdom: object | None,
    phylum: object | None = None,
    lineage_class: object | None = None,
    *,
    rules: OrganismClassificationRules | None = None,
) -> str | None:
    """Return the cellularity label for a single taxonomy triple.

    Parameters
    ----------
    superkingdom, phylum, lineage_class:
        Raw taxonomy values.  They are normalised using :func:`normalize`
        before matching against the configured rule buckets.
    rules:
        Classification rules.  Defaults to :data:`DEFAULT_RULES`.

    Returns
    -------
    str | None
        Matching cellularity label or ``None`` when the lineage does not
        provide enough information to classify the organism.
    """

    rules = rules or DEFAULT_RULES

    sk = normalize(superkingdom)
    ph = normalize(phylum)
    cl = normalize(lineage_class)

    if sk in rules.virus_superkingdoms or ph in rules.virus_phyla:
        return rules.label_viruses

    if cl and cl in rules.multicellular_classes:
        return rules.label_multicellular
    if cl and cl in rules.unicellular_classes:
        return rules.label_unicellular

    if ph and ph in rules.multicellular_phyla:
        return rules.label_multicellular
    if ph and ph in rules.unicellular_phyla:
        return rules.label_unicellular

    if sk and sk in rules.unicellular_superkingdoms:
        return rules.label_unicellular

    return None


def classify_record(
    record: Mapping[str, object] | pd.Series,
    *,
    rules: OrganismClassificationRules | None = None,
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str = "lineage_class",
) -> str | None:
    """Classify a mapping containing lineage information."""

    rules = rules or DEFAULT_RULES

    getter: Mapping[str, object]
    if isinstance(record, Mapping):
        getter = record
    elif hasattr(record, "get"):
        getter = record  # pandas Series implements ``get``
    else:
        raise TypeError("record must be a mapping or pandas Series")

    superkingdom = getter.get(superkingdom_col)
    phylum = getter.get(phylum_col)
    lineage_class = getter.get(class_col)
    return classify_by_lineage(
        superkingdom,
        phylum,
        lineage_class,
        rules=rules,
    )


def _extract_column(df: pd.DataFrame, column: str | None) -> list[object | None]:
    if column is None:
        return [None] * len(df)
    if column not in df.columns:
        raise KeyError(f"Column '{column}' not found in dataframe")
    return df[column].tolist()


def _classify_dataframe(
    df: pd.DataFrame,
    *,
    rules: OrganismClassificationRules,
    superkingdom_col: str | None,
    phylum_col: str | None,
    class_col: str | None,
) -> pd.Series:
    supers = _extract_column(df, superkingdom_col)
    phyla = _extract_column(df, phylum_col)
    classes = _extract_column(df, class_col)
    values = [
        classify_by_lineage(s, p, c, rules=rules)
        for s, p, c in zip(supers, phyla, classes, strict=True)
    ]
    return pd.Series(pd.array(values, dtype="string"), index=df.index)


def add_cellularity(
    df: pd.DataFrame,
    *,
    rules: OrganismClassificationRules | None = None,
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str | None = "lineage_class",
    target_col: str = "cellularity",
) -> pd.DataFrame:
    """Append a cellularity column derived from taxonomy lineage fields."""

    rules = rules or DEFAULT_RULES
    result = df.copy()
    result[target_col] = _classify_dataframe(
        result,
        rules=rules,
        superkingdom_col=superkingdom_col,
        phylum_col=phylum_col,
        class_col=class_col,
    )
    return result


def add_cellularity_smart(
    df: pd.DataFrame,
    *,
    rules: OrganismClassificationRules | None = None,
    superkingdom_col: str = "lineage_superkingdom",
    phylum_col: str = "lineage_phylum",
    class_col: str | None = "lineage_class",
    target_col: str = "organism_type",
) -> pd.DataFrame:
    """Fill ``target_col`` with cellularity, preserving existing labels."""

    rules = rules or DEFAULT_RULES
    result = df.copy()
    derived = _classify_dataframe(
        result,
        rules=rules,
        superkingdom_col=superkingdom_col,
        phylum_col=phylum_col,
        class_col=class_col,
    )

    if target_col in result.columns:
        existing = result[target_col].astype("string")
        keep_mask = existing.map(normalize).notna()
        result[target_col] = existing.where(keep_mask, derived)
    else:
        result[target_col] = derived
    return result


__all__ = [
    "OrganismClassificationRules",
    "add_cellularity",
    "add_cellularity_smart",
    "classify_by_lineage",
    "classify_record",
    "normalize",
]
