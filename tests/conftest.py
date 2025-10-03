from __future__ import annotations

import csv
import datetime as dt
import os
import random
import time
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


FROZEN_UTC = datetime(2020, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
FROZEN_TIMESTAMP = FROZEN_UTC.timestamp()
FROZEN_NAIVE = FROZEN_UTC.replace(tzinfo=None)


def _fix_seed(
    seed: int = 42, *, monkeypatch: pytest.MonkeyPatch | None = None
) -> None:
    if monkeypatch is None:
        os.environ["PYTHONHASHSEED"] = str(seed)
    else:
        monkeypatch.setenv("PYTHONHASHSEED", str(seed))
    random.seed(seed)
    np.random.seed(seed)


@pytest.fixture(autouse=True)
def deterministic_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Provide deterministic environment defaults for every test."""

    _fix_seed(monkeypatch=monkeypatch)
    monkeypatch.setenv("TZ", "UTC")
    monkeypatch.setenv("HOME", str(tmp_path))

    class FrozenDateTime(dt.datetime):
        @classmethod
        def now(cls, tz: dt.tzinfo | None = None) -> dt.datetime:
            if tz is None:
                return FROZEN_NAIVE
            return FROZEN_UTC.astimezone(tz)

        @classmethod
        def utcnow(cls) -> dt.datetime:
            return FROZEN_NAIVE

        @classmethod
        def today(cls) -> dt.datetime:
            return cls.now()

    monkeypatch.setattr(time, "time", lambda: FROZEN_TIMESTAMP)
    monkeypatch.setattr(dt, "datetime", FrozenDateTime)
    monkeypatch.setattr("datetime.datetime", FrozenDateTime)


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
    return FROZEN_UTC.isoformat()
