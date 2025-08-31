"""Utility functions for reading input data.

This module provides helpers to read CSV/TSV files with automatic
separator and encoding detection. Column names are normalised to lower
case and mapped to canonical aliases used by the classifier.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Generator, Iterable, Optional
import csv
import logging

import pandas as pd

# ---------------------------------------------------------------------------
# Column aliases
# ---------------------------------------------------------------------------

PT_ALIASES: Dict[str, Iterable[str]] = {
    "pubmed_pt": ["pubmed.publicationtype", "pubmed_publicationtype", "pubmed_pt"],
    "openalex_pt": [
        "openalex.publicationtypes",
        "openalex_publicationtypes",
        "openalex_pt",
    ],
    "scholar_pt": ["scholar.publicationtypes", "scholar_publicationtypes", "scholar_pt"],
    "crossref_type": ["crossref.type", "crossref_type"],
}

MESH_ALIASES: Dict[str, Iterable[str]] = {
    "pubmed_mesh_desc": [
        "pubmed.mesh_descriptors",
        "pubmed_mesh_descriptors",
        "pubmed_mesh_desc",
    ],
    "pubmed_mesh_qual": [
        "pubmed.mesh_qualifiers",
        "pubmed_mesh_qualifiers",
        "pubmed_mesh_qual",
    ],
    "openalex_mesh_desc": [
        "openalex.meshdescriptors",
        "openalex_meshdescriptors",
        "openalex_mesh_desc",
    ],
    "openalex_mesh_qual": [
        "openalex.meshqualifiers",
        "openalex_meshqualifiers",
        "openalex_mesh_qual",
    ],
}

ID_COLUMNS = [
    "document_chembl_id",
    "doi",
    "pubmed_id",
    "openalex.paperid",
]

ALL_ALIASES: Dict[str, Iterable[str]] = {**PT_ALIASES, **MESH_ALIASES}

# ---------------------------------------------------------------------------
# Reading helpers
# ---------------------------------------------------------------------------


def _detect_encoding(path: str) -> str:
    """Attempt to detect file encoding.

    Tries UTF-8 first and falls back to latin-1 on failure. This keeps the
    dependency footprint minimal while handling the majority of files.
    """
    try:
        with open(path, "r", encoding="utf-8") as fh:
            fh.read(1024)
        return "utf-8"
    except UnicodeDecodeError:
        logging.warning("Failed to decode %s as UTF-8, falling back to latin-1", path)
        return "latin-1"


def read_chunks(
    path: str,
    chunk_size: int = 1000,
    sep: Optional[str] = None,
) -> Generator[pd.DataFrame, None, None]:
    """Yield chunks of the input CSV/TSV file.

    Parameters
    ----------
    path:
        Path to the input CSV/TSV file.
    chunk_size:
        Number of rows per chunk.
    sep:
        Optional column separator. If ``None`` the separator will be inferred
        using the Python engine of :func:`pandas.read_csv`.
    """

    encoding = _detect_encoding(path)
    logging.info("Reading %s with encoding %s", path, encoding)

    reader = pd.read_csv(
        path,
        sep=sep,
        engine="python" if sep is None else None,
        chunksize=chunk_size,
        dtype=str,
        encoding=encoding,
        keep_default_na=False,
    )

    for chunk in reader:
        # Normalise column names
        chunk.columns = [c.lower().strip() for c in chunk.columns]
        for canonical, aliases in ALL_ALIASES.items():
            for alias in aliases:
                if alias.lower() in chunk.columns and canonical not in chunk.columns:
                    chunk.rename(columns={alias.lower(): canonical}, inplace=True)
        yield chunk


__all__ = ["read_chunks", "ID_COLUMNS", "PT_ALIASES", "MESH_ALIASES"]
