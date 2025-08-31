"""Text normalisation utilities.

This module performs normalisation of PublicationType and MeSH strings.
The main entry point :func:`split_and_canon` converts raw strings into a
list of canonical tokens:

* lower-case transformation
* replacement of known synonyms to a canonical form
* tokenisation by common delimiters
* trimming of whitespace and removal of duplicates while preserving
  order
"""
from __future__ import annotations

import re
from typing import Dict, List

# Regular expression matching standard delimiters: pipe, semicolon, comma and
# slash. Spaces are handled by stripping and thus remain part of the token to
# preserve multi-word phrases such as ``randomized controlled trial``.
_DELIMS_RE = re.compile(r"[|;,/]+")

# ---------------------------------------------------------------------------
# Canonicalisation tables
# ---------------------------------------------------------------------------

# Mapping of various strings to canonical forms. The keys are matched in the
# lower-cased input string before tokenisation.
CANON_TABLE: Dict[str, str] = {
    # Review synonyms
    "review article": "review",
    "review-article": "review",
    "mini review": "review",
    "mini-review": "review",
    "literature review": "review",
    "umbrella review": "review",
    "narrative review": "review",
    "state-of-the-art review": "review",
    "integrative review": "review",
    "critical review": "review",
    "brief review": "review",
    "comprehensive review": "review",
    "evidence synthesis": "review",
    "qualitative review": "review",
    "overview of reviews": "review",
    "scoping review": "systematic review",
    "systematic literature review": "systematic review",
    "slr": "systematic review",
    "meta analysis": "meta-analysis",
    "meta-analysis": "meta-analysis",
    "network meta-analysis": "meta-analysis",
    "meta-synthesis": "meta-analysis",
    "journalarticle": "journal-article",
    "journal article": "journal-article",
    "journal-article": "journal-article",
    "research article": "journal-article",
    "original article": "journal-article",
    "brief communication": "journal-article",
    "rapid communication": "journal-article",
}

# Additional replacements can be added by users of the library by extending
# this dictionary.


def split_and_canon(value: str | None, canon_table: Dict[str, str] | None = None) -> List[str]:
    """Split a string into normalised canonical tokens.

    Parameters
    ----------
    value:
        Raw string value from the CSV file.
    canon_table:
        Optional synonym mapping table. If omitted ``CANON_TABLE`` is used.

    Returns
    -------
    list[str]
        List of unique canonical tokens in the order of appearance.
    """

    if value is None or value == "" or (isinstance(value, float) and value != value):
        return []

    canon_table = canon_table or CANON_TABLE

    # Lower-case
    text = str(value).lower()

    # Replace known synonyms before splitting. Word boundaries ensure that
    # only complete tokens are substituted.
    for pattern, replacement in canon_table.items():
        text = re.sub(rf"\b{re.escape(pattern)}\b", replacement, text)

    # Tokenise
    tokens = [t.strip() for t in _DELIMS_RE.split(text) if t.strip()]

    # Deduplicate while preserving order
    seen = set()
    result: List[str] = []
    for token in tokens:
        if token not in seen:
            seen.add(token)
            result.append(token)
    return result


__all__ = ["split_and_canon", "CANON_TABLE"]
