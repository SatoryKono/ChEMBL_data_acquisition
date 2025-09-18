"""Test configuration fixtures.

This module provides common pytest fixtures used across the test suite.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.config import ApiCfg, Config


@pytest.fixture()
def cfg() -> Config:
    """Return a baseline :class:`~library.config.Config` instance for tests.

    The configuration requires a valid ``api.user_agent`` value; the fixture
    supplies a deterministic placeholder address so that individual tests do
    not need to construct :class:`Config` instances manually.
    """

    return Config(api=ApiCfg(user_agent="test@example.org"))


@pytest.fixture()
def duplicate_document_ids() -> list[str]:
    """Return sample document IDs including duplicates for testing."""

    return ["CHEMBL1", "CHEMBL1", "CHEMBL2"]


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture()
def activities_sample_df() -> pd.DataFrame:
    """Load the minimal activities CSV as a :class:`~pandas.DataFrame`."""

    return pd.read_csv(
        DATA_DIR / "activities_sample.csv",
        dtype={
            "activity_id": str,
            "molecule_chembl_id": str,
            "assay_chembl_id": str,
            "standard_type": str,
        },
    )


@pytest.fixture()
def documents_sample_df() -> pd.DataFrame:
    """Load the minimal documents CSV as a :class:`~pandas.DataFrame`."""

    return pd.read_csv(DATA_DIR / "documents_sample.csv")
