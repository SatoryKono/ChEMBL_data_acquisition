"""Normalization utilities for publication data."""
from __future__ import annotations

import re
from typing import Dict, List, Mapping, Sequence, Tuple

from .constants import PUBLICATION_TYPE_MAPPING


def normalize_publication_types(
    raw: Mapping[str, Sequence[str]],
) -> Dict[str, List[str]]:
    """Normalize publication type labels from multiple sources.

    Parameters
    ----------
    raw:
        Mapping of source name to list of raw publication type strings.

    Returns
    -------
    dict
        Mapping of source name to list of normalized publication type strings.
    """
    normalized: Dict[str, List[str]] = {}
    for source, pts in raw.items():
        mapping = PUBLICATION_TYPE_MAPPING.get(source, {})
        norm_pts = []
        for pt in pts:
            mapped = mapping.get(pt)
            if mapped:
                norm_pts.append(mapped)
        if norm_pts:
            normalized[source] = norm_pts
    return normalized


def _clean(text: str) -> str:
    """Normalize text by lowering case and stripping punctuation."""
    lowered = text.lower()
    cleaned = re.sub(r"[\-\._]", " ", lowered)
    cleaned = re.sub(r"[^a-z0-9\s]", "", cleaned)
    return cleaned.strip()


def normalize_mesh(
    descriptors: Mapping[str, Sequence[str]],
    qualifiers: Mapping[str, Sequence[str]],
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Normalize MeSH descriptors and qualifiers.

    Parameters
    ----------
    descriptors, qualifiers:
        Mapping of source name to list of descriptor/qualifier strings.

    Returns
    -------
    tuple(dict, dict)
        Two dictionaries with normalized descriptors and qualifiers.
    """
    norm_desc: Dict[str, List[str]] = {}
    norm_qual: Dict[str, List[str]] = {}
    for source, items in descriptors.items():
        norm_desc[source] = [_clean(d) for d in items]
    for source, items in qualifiers.items():
        norm_qual[source] = [_clean(q) for q in items]
    return norm_desc, norm_qual
