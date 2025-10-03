from __future__ import annotations

import csv
import json
import os
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import pytest

from library.config import Config


def _fix_seed(seed: int = 42) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


@pytest.fixture(autouse=True)
def deterministic_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Provide deterministic environment defaults for every test."""

    _fix_seed()
    monkeypatch.setenv("TZ", "UTC")
    monkeypatch.setenv("HOME", str(tmp_path))


@pytest.fixture(autouse=True)
def disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Disallow outbound HTTP requests during tests."""

    try:
        import requests
    except ModuleNotFoundError:  # pragma: no cover - requests optional in env
        return

    def deny(self, method, url, *args, **kwargs):  # type: ignore[override]
        raise AssertionError("External network access is disabled during tests")

    monkeypatch.setattr("requests.sessions.Session.request", deny)


@pytest.fixture()
def cfg() -> Config:
    """Return a baseline :class:`~library.config.Config` instance for tests."""

    config = Config()
    config.api.user_agent = "test@example.org"
    config.sources.pubchem.user_agent = "pubchem-contact@example.org"
    return config


@pytest.fixture()
def sample_input_csv(tmp_path: Path) -> Path:
    path = tmp_path / "input.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["molecule_chembl_id", "name", "smiles"])
        writer.writerow(["CHEMBL1", "Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"])
        writer.writerow(["CHEMBL2", "Unknown", ""])
    return path


@pytest.fixture()
def snapshot_resource() -> Path:
    return Path(__file__).parent / "resources"


@pytest.fixture()
def make_dataframe() -> Iterable[pd.DataFrame]:
    """Return a factory creating deterministic test dataframes."""

    def _factory(rows: Iterable[dict[str, object]]) -> pd.DataFrame:
        return pd.DataFrame(rows)

    return _factory


@pytest.fixture()
def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
