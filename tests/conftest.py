from __future__ import annotations

import csv
import datetime as dt
import os
import random
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Sequence

import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
import pytest

from library.config import Config
from library.orchestration import ETLContext


FROZEN_UTC = datetime(2020, 1, 1, 0, 0, 0, tzinfo=timezone.utc)
FROZEN_TIMESTAMP = FROZEN_UTC.timestamp()
FROZEN_NAIVE = FROZEN_UTC.replace(tzinfo=None)


def _fix_seed(
    seed: int = 42, *, monkeypatch: pytest.MonkeyPatch | None = None
) -> None:
    """Reset the pseudo-random generators to a deterministic state."""

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

    class FrozenDate(dt.date):
        @classmethod
        def today(cls) -> dt.date:
            return FROZEN_NAIVE.date()

    monkeypatch.setattr(time, "time", lambda: FROZEN_TIMESTAMP)
    monkeypatch.setattr(time, "time_ns", lambda: int(FROZEN_TIMESTAMP * 1_000_000_000))
    monkeypatch.setattr(dt, "datetime", FrozenDateTime)
    monkeypatch.setattr("datetime.datetime", FrozenDateTime)
    monkeypatch.setattr(dt, "date", FrozenDate)
    monkeypatch.setattr("datetime.date", FrozenDate)


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
def make_fallback_doi_csv(tmp_path: Path) -> Callable[..., Path]:
    """Return a factory producing deterministic fallback DOI CSV files."""

    def _factory(
        rows: Iterable[tuple[str, str]],
        *,
        columns: Sequence[str] = ("PMID", "DOI"),
        filename: str = "fallback_doi.csv",
    ) -> Path:
        path = tmp_path / filename
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(list(columns))
            for pmid, doi in rows:
                writer.writerow([pmid, doi])
        return path

    return _factory


@pytest.fixture()
def fallback_doi_csv(
    make_fallback_doi_csv: Callable[..., Path]
) -> Path:
    """Provide a default fallback DOI CSV used across tests."""

    default_rows = (
        ("123456", "10.1000/default-one"),
        ("789012", "10.1000/default-two"),
    )
    return make_fallback_doi_csv(default_rows)


@pytest.fixture()
def make_dataframe() -> Iterable[pd.DataFrame]:
    """Return a factory creating deterministic test dataframes."""

    def _factory(rows: Iterable[dict[str, object]]) -> pd.DataFrame:
        return pd.DataFrame(rows)

    return _factory


@pytest.fixture()
def utc_now_iso() -> str:
    return FROZEN_UTC.isoformat()


class _StubChemblClient:
    def __init__(self, registry: dict[str, list["_StubChemblClient"]], *args, **kwargs):
        self._registry = registry
        self._registry.setdefault("created", []).append(self)
        self._registry.setdefault("closed", [])
        self._closed = False

    def close(self) -> None:
        if not self._closed:
            self._registry.setdefault("closed", []).append(self)
            self._closed = True

    def __enter__(self) -> "_StubChemblClient":  # pragma: no cover - trivial helper
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:  # pragma: no cover - trivial helper
        self.close()
        return False


@pytest.fixture()
def stub_etl_context() -> "StubETLFactory":
    registry: dict[str, list[_StubChemblClient]] = {"created": [], "closed": []}

    class StubETLFactory:
        def __call__(self, cfg: Config, **kwargs: object) -> ETLContext:
            return ETLContext(
                cfg,
                chembl_client_factory=lambda *args, **kw: _StubChemblClient(
                    registry, *args, **kw
                ),
                **kwargs,
            )

        @property
        def created_clients(self) -> list[_StubChemblClient]:
            return registry["created"]

        @property
        def closed_clients(self) -> list[_StubChemblClient]:
            return registry["closed"]

    return StubETLFactory()
