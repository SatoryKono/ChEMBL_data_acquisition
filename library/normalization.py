"""Simple DataFrame normalisation helpers."""

from __future__ import annotations

import pandas as pd


def normalize_activities(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise ChEMBL activity records.

    Parameters
    ----------
    df:
        Raw activity DataFrame.

    Returns
    -------
    pandas.DataFrame
        DataFrame with standardised dtypes.
    """
    return df.convert_dtypes()


def normalize_assays(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise assay records returned by the ChEMBL API."""
    return df.convert_dtypes()


def normalize_documents(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise document metadata."""
    return df.convert_dtypes()


def normalize_targets(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise target information."""
    return df.convert_dtypes()


def normalize_testitems(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise test item (compound) information."""
    return df.convert_dtypes()
