"""Term dictionaries and normalization utilities for publication classification.

This module defines term sets used for classification and provides functions to
normalise publication type fields.
"""
from __future__ import annotations

import re
from typing import Iterable, List

# Term dictionaries -----------------------------------------------------------------

REVIEW_TERMS: set[str] = {
    "review",
    "book review",
    "systematic review",
    "meta-analysis",
    "scoping review",
    "umbrella review",
    "mini-review",
    "correction",
    "survey",
    "overview",
    "state of the art",
    "bibliometric analysis",
    "editorial",
    "corrected and republished article",
    "address",
    "lecture",
    "news",
    "historical article",
    "comment",
    "video-audio media",
    "retracted publication",
    # OpenAlex specific
    "review-article",

}

EXPERIMENTAL_TERMS: set[str] = {
    "comparative study",
    "evaluation study",
    "validation study",
    "case-control study",
    "cohort study",
    "proceedings-article",
    "posted-content",
    "preprint",
}

UNKNOWN_TERMS: set[str] = {
    "research support, non-u.s. gov't",
    "research support, n.i.h., extramural",
    "research support, u.s. gov't, non-p.h.s.",
    "research support, u.s. gov't, p.h.s.",
    "research support, n.i.h., intramural",
    "clinical trial",
    "randomized controlled trial",
    "clinical study",
    "case-control study",
    "cohort study",
    "note",
    "erratum",
    "data paper",
    "perspective",
    "opinion",
    "short survey",
}

# Normalisation utilities ------------------------------------------------------------

# Map various synonyms and fused forms to canonical tokens.
_SYNONYMS: dict[str, str] = {
    "review article": "review",
    "review-article": "review",
    "mini review": "review",
    "mini-review": "review",
    "meta analysis": "meta-analysis",
    "meta-analysis": "meta-analysis",
    "journal article": "journal-article",
    "journal-article": "journal-article",
    "journalarticle": "journal-article",
    "slr": "systematic review",
    "state-of-the-art": "state of the art",
}

_DELIMITERS_RE = re.compile(r"[|;,/]")


def _normalise_token(token: str) -> str:
    """Return canonical representation of a token.

    Parameters
    ----------
    token:
        Raw token extracted from a publication type field.

    Returns
    -------
    str
        Canonical form or an empty string if token is empty after processing.
    """
    token = token.strip().lower()
    if not token:
        return ""
    token = _SYNONYMS.get(token, token)
    return token


def parse_terms(value: object) -> List[str]:
    """Split and normalise a publication type field.

    Parameters
    ----------
    value:
        Field value which may be a string or missing (NaN/None).

    Returns
    -------
    list[str]
        Sorted list of unique canonical tokens. Empty list if no valid tokens.
    """
    if not isinstance(value, str):
        return []

    parts = _DELIMITERS_RE.split(value)
    tokens = [_normalise_token(p) for p in parts]
    tokens = [t for t in tokens if t]

    # Remove duplicates while preserving order, then sort for deterministic output
    unique = []
    seen = set()
    for t in tokens:
        if t not in seen:
            unique.append(t)
            seen.add(t)

    return sorted(unique)
