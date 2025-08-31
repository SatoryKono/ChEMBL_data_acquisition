from __future__ import annotations

"""Normalization helpers for publication types and MeSH terms."""

from collections import OrderedDict
import re
from typing import Callable, List

import pandas as pd

# Aliases that map to ``review``
REVIEW_ALIASES: dict[str, str] = {
    "review": "review",
    "systematic review": "review",
    "meta-analysis": "review",
    "meta analysis": "review",
    "literature review": "review",
    "mini review": "review",
    "narrative review": "review",
    "scoping review": "review",
    "umbrella review": "review",
    "rapid review": "review",
    "evidence synthesis": "review",
}


def _split_tokens(value: str, separators: str) -> List[str]:
    pattern = f"[{re.escape(separators)}]"
    tokens = re.split(pattern, value)
    cleaned = [tok.strip() for tok in tokens if tok and tok.strip()]
    # deduplicate preserving order
    return list(OrderedDict.fromkeys(cleaned))


def split_and_normalize(value: object, separators: str, *, canonicalizer: Callable[[str], str] | None = None) -> List[str]:
    """Split string ``value`` by ``separators`` and normalize.

    Parameters
    ----------
    value:
        Raw string value.
    separators:
        String with separator characters.
    canonicalizer:
        Optional callable applied to each token.

    Returns
    -------
    list[str]
        Normalized, deduplicated tokens.
    """
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    tokens = _split_tokens(str(value).lower(), separators)
    if canonicalizer:
        tokens = [canonicalizer(tok) for tok in tokens]
    return tokens


def is_review(types_set: set[str]) -> bool:
    """Return ``True`` if any token denotes a review article."""
    normalized = {REVIEW_ALIASES.get(tok, tok) for tok in types_set}
    return "review" in normalized


def normalize_publication_types(df: pd.DataFrame, separators: str) -> pd.DataFrame:
    """Add normalized publication type columns to ``df``."""
    columns = {
        "PubMed.PublicationType": "norm.PubMed.PublicationType",
        "OpenAlex.PublicationTypes": "norm.OpenAlex.PublicationTypes",
        "scholar.PublicationTypes": "norm.scholar.PublicationTypes",
    }
    for src, dst in columns.items():
        df[dst] = df[src].apply(lambda x: split_and_normalize(x, separators))
    return df


def canonicalize_mesh(term: str) -> str:
    """Canonically normalise a MeSH term.

    This is a heuristic implementation intended for matching to the
    ``experimental_mesh`` table.
    """
    term = term.lower().replace("&", "and")
    term = re.sub(r"[–‐−]", "-", term)  # unify dashes
    term = re.sub(r"\(.*?\)", "", term)  # remove bracketed noise
    term = term.replace(".", "")
    term = term.strip()
    # naive singularisation
    if term.endswith("ies"):
        term = term[:-3] + "y"
    elif term.endswith("ses"):
        term = term[:-2]
    elif term.endswith("s") and not term.endswith("ss"):
        term = term[:-1]
    return term


def normalize_mesh_fields(df: pd.DataFrame, separators: str) -> pd.DataFrame:
    """Add normalised MeSH columns to ``df``."""
    columns = {
        "PubMed.MeSH_Descriptors": "norm.PubMed.MeSH_Descriptors",
        "PubMed.MeSH_Qualifiers": "norm.PubMed.MeSH_Qualifiers",
        "OpenAlex.MeshDescriptors": "norm.OpenAlex.MeshDescriptors",
        "OpenAlex.MeshQualifiers": "norm.OpenAlex.MeshQualifiers",
    }
    for src, dst in columns.items():
        df[dst] = df[src].apply(lambda x: split_and_normalize(x, separators, canonicalizer=canonicalize_mesh))
    return df


def collect_mesh_terms(row: pd.Series) -> list[str]:
    """Collect unique MeSH terms from normalised columns in ``row``."""
    fields = [
        "norm.PubMed.MeSH_Descriptors",
        "norm.PubMed.MeSH_Qualifiers",
        "norm.OpenAlex.MeshDescriptors",
        "norm.OpenAlex.MeshQualifiers",
    ]
    all_terms: list[str] = []
    for f in fields:
        all_terms.extend(row.get(f, []))
    return list(OrderedDict.fromkeys(all_terms))
