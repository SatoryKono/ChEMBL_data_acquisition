"""Wrappers around :mod:`schemas.normalize` helpers."""

from __future__ import annotations

import pandas as pd
from schemas.normalize import (
    normalize_activities as _normalize_activities,
    normalize_assays as _normalize_assays,
    normalize_documents as _normalize_documents,
    normalize_targets as _normalize_targets,
    normalize_testitems as _normalize_testitems,
)


def normalize_activities(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise ChEMBL activity records.

    This is a thin wrapper around :func:`schemas.normalize.normalize_activities`.

    Parameters
    ----------
    df:
        Raw activity DataFrame.

    Returns
    -------
    pandas.DataFrame
        Normalised copy of ``df``.
    """

    return _normalize_activities(df)


def normalize_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise assay records returned by the ChEMBL API.

    This delegates to :func:`schemas.normalize.normalize_assays`.
    """

    return _normalize_assays(df)


def normalize_documents(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise document metadata.

    This delegates to :func:`schemas.normalize.normalize_documents`.
    """

    return _normalize_documents(df)


def normalize_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise target information.

    This delegates to :func:`schemas.normalize.normalize_targets`.
    """

    return _normalize_targets(df)


def normalize_testitems(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise test item (compound) information.

    This delegates to :func:`schemas.normalize.normalize_testitems`.
    """

    return _normalize_testitems(df)


__all__ = [
    "normalize_activities",
    "normalize_assays",
    "normalize_documents",
    "normalize_targets",
    "normalize_testitems",
]
