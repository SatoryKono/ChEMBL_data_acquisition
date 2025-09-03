from __future__ import annotations

"""Utilities for working with experimental MeSH probabilities."""

from pathlib import Path

from .io import read_mesh_probabilities
from .normalize import canonicalize_mesh

MeshProbMap = dict[str, float]


__all__ = ["MeshProbMap", "load_mesh_probabilities"]


def load_mesh_probabilities(path: Path) -> MeshProbMap:
    """Load experimental MeSH probabilities from ``path``."""
    raw = read_mesh_probabilities(path)
    return {canonicalize_mesh(term): prob for term, prob in raw.items()}
