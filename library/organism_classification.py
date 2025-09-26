"""Utilities for deriving organism cellularity information from taxonomy."""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

DEFAULT_TYPE = "-"
VIRAL_SUPERKINGDOMS = {"viruses"}
UNICELLULAR_SUPERKINGDOMS = {"bacteria", "archaea"}
VIRAL_KEYWORDS: tuple[str, ...] = ("virus", "phage")
UNICELLULAR_PHYLUMS = {
    "amoebozoa",
    "bacillati",
    "cryptista",
    "discoba",
    "fungi",
    "haptista",
    "metamonada",
    "methanobacteriati",
    "myzozoa",
    "pseudomonadati",
    "sar",
    "thermotogati",
}


def _normalise_token(value: object) -> str:
    if isinstance(value, str):
        return value.strip()
    if value is None:
        return ""
    if pd.isna(value):  # type: ignore[arg-type]
        return ""
    return str(value).strip()


def _is_viral(genus: str, superkingdom: str) -> bool:
    if superkingdom.lower() in VIRAL_SUPERKINGDOMS:
        return True
    genus_lower = genus.lower()
    return any(token in genus_lower for token in VIRAL_KEYWORDS)


def _is_unicellular(superkingdom: str, phylum: str, lineage_class: str) -> bool:
    if superkingdom.lower() in UNICELLULAR_SUPERKINGDOMS:
        return True
    tokens: Iterable[str] = (phylum.lower(), lineage_class.lower())
    return any(token in UNICELLULAR_PHYLUMS for token in tokens if token)


def add_cellularity_smart(
    df: pd.DataFrame,
    *,
    genus_col: str = "genus",
    superkingdom_col: str = "superkingdom",
    phylum_col: str = "phylum",
    lineage_class_col: str = "lineage_class",
    taxon_id_col: str = "taxon_id",
) -> pd.DataFrame:
    """Return organism cellularity derived from taxonomy hints.

    Parameters
    ----------
    df:
        DataFrame containing organism taxonomy columns.
    genus_col, superkingdom_col, phylum_col, lineage_class_col, taxon_id_col:
        Column names providing genus and taxonomy attributes. Missing columns
        are treated as empty strings.

    Returns
    -------
    pandas.DataFrame
        Two-column frame with ``genus`` and the inferred ``type`` suitable for
        :func:`library.target_postprocessing.finalise_targets`.
    """

    if genus_col not in df.columns:
        return pd.DataFrame({"genus": pd.Series(dtype="string"), "type": pd.Series(dtype="string")})

    genus = df[genus_col].astype("string")
    superkingdom = df.get(superkingdom_col, pd.Series(index=df.index, dtype="string")).astype("string")
    phylum = df.get(phylum_col, pd.Series(index=df.index, dtype="string")).astype("string")
    lineage_class = df.get(lineage_class_col, pd.Series(index=df.index, dtype="string")).astype("string")

    records = []
    for g, sk, ph, lc in zip(genus, superkingdom, phylum, lineage_class, strict=False):
        g_norm = _normalise_token(g)
        if not g_norm:
            continue
        sk_norm = _normalise_token(sk)
        ph_norm = _normalise_token(ph)
        lc_norm = _normalise_token(lc)
        if _is_viral(g_norm, sk_norm):
            inferred = "Viruses"
        elif _is_unicellular(sk_norm, ph_norm, lc_norm):
            inferred = "Unicellular organism"
        elif sk_norm:
            inferred = "Multicellular organism"
        else:
            inferred = DEFAULT_TYPE
        records.append((g_norm, inferred))

    if not records:
        return pd.DataFrame({"genus": pd.Series(dtype="string"), "type": pd.Series(dtype="string")})

    priority = {"Viruses": 3, "Unicellular organism": 2, "Multicellular organism": 1, DEFAULT_TYPE: 0}
    df_records = pd.DataFrame(records, columns=["genus", "type"])
    df_records = df_records.sort_values(by="type", key=lambda s: s.map(priority).fillna(0), ascending=False)
    deduplicated = df_records.drop_duplicates(subset="genus", keep="first")
    deduplicated["type"] = deduplicated["type"].astype("string")
    return deduplicated.reset_index(drop=True)

