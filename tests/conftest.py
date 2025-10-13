from __future__ import annotations

import csv
import datetime as dt
import os
import random
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from datetime import datetime
from importlib import import_module
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# ruff: noqa: E402  # test utilities modify sys.path before project imports
import numpy as np
import pandas as pd
import pytest
import yaml

try:  # pragma: no cover - optional dependency shim
    import pandera.pandas as pa  # type: ignore  # noqa: F401
except ModuleNotFoundError:  # pragma: no cover - fallback for Pandera >=0.20
    import types

    pandas_model = import_module("pandera.api.pandas.model")
    shim = types.ModuleType("pandera.pandas")
    for attr in dir(pandas_model):
        if attr.startswith("__"):
            continue
        setattr(shim, attr, getattr(pandas_model, attr))
    sys.modules["pandera.pandas"] = shim

from config.paths import DICTIONARY_DIR
from library.config import Config
from library.orchestration import ETLContext
from library.resources import dictionaries as dictionary_resources
from tests.helpers.reporting_policy import ReportingPolicyPlugin

FROZEN_UTC = datetime(2020, 1, 1, 0, 0, 0, tzinfo=dt.UTC)
FROZEN_TIMESTAMP = FROZEN_UTC.timestamp()
FROZEN_NAIVE = FROZEN_UTC.replace(tzinfo=None)


def _fix_seed(seed: int = 42, *, monkeypatch: pytest.MonkeyPatch | None = None) -> None:
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
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    monkeypatch.setenv("TMP", str(tmp_path))
    monkeypatch.setenv("TEMP", str(tmp_path))

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


def _is_truthy(value: str) -> bool:
    return value.lower() in {"1", "true", "yes", "on"}


def _allow_dictionary_relaxation() -> bool:
    flag = os.getenv("CHEMBL_DA_ALLOW_DICT_RELAX")
    if flag is None:
        return False
    return _is_truthy(flag.strip())


@pytest.fixture(scope="session", autouse=True)
def relax_dictionary_manifest_checks() -> None:
    """Allow dictionary metadata lookup even when bundled samples diverge."""

    try:
        dictionary_resources.list_resources()
    except dictionary_resources.DictionaryManifestError:
        manifest_path = DICTIONARY_DIR / "manifest.yaml"
        try:
            manifest_data = (
                yaml.safe_load(manifest_path.read_text(encoding="utf-8")) or {}
            )
        except OSError:
            yield
            return

        resources = manifest_data.get("resources")
        if not isinstance(resources, dict):
            yield
            return

        original_get_resource = dictionary_resources.get_resource

        def _safe_get_resource(name: str, *, base_dir: Path | None = None):
            try:
                return original_get_resource(name, base_dir=base_dir)
            except dictionary_resources.DictionaryManifestError:
                if not _allow_dictionary_relaxation():
                    raise
                entry = resources.get(name)
                if not isinstance(entry, dict):
                    raise
                path_value = entry.get("path")
                version = entry.get("version")
                sha256 = entry.get("sha256")
                generator = entry.get("generator", "")
                if (
                    not isinstance(path_value, str)
                    or not isinstance(version, str)
                    or not isinstance(sha256, str)
                ):
                    raise
                root = Path(base_dir) if base_dir is not None else DICTIONARY_DIR
                resolved_path = (root / path_value).resolve()
                return dictionary_resources.DictionaryResource(
                    name=name,
                    relative_path=Path(path_value),
                    path=resolved_path,
                    version=version,
                    sha256=sha256,
                    generator=Path(generator),
                )

        dictionary_resources.get_resource = _safe_get_resource
        try:
            yield
        finally:
            dictionary_resources.get_resource = original_get_resource
    else:
        yield


@pytest.fixture(autouse=True)
def disable_network(monkeypatch: pytest.MonkeyPatch) -> None:
    """Disallow outbound HTTP requests during tests."""

    try:
        import requests  # noqa: F401
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
def stub_etl_context(cfg: Config):
    """Return an :class:`ETLContext` wired with stubbed Chembl clients."""

    class _StubClient:
        def __init__(self) -> None:
            self.close_calls = 0

        def close(self) -> None:
            self.close_calls += 1

    created_clients: list[_StubClient] = []

    def _factory(_context: ETLContext) -> _StubClient:
        client = _StubClient()
        created_clients.append(client)
        return client

    context = ETLContext(cfg, chembl_client_factory=_factory)
    try:
        yield context, created_clients
    finally:
        context.close()


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
def fallback_doi_csv(make_fallback_doi_csv: Callable[..., Path]) -> Path:
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


_REPORTING_POLICY_PLUGIN_ATTR = "_reporting_policy_plugin"


def pytest_configure(config: pytest.Config) -> None:
    """Register the reporting policy plugin and prepare report directories."""

    plugin = getattr(config, _REPORTING_POLICY_PLUGIN_ATTR, None)
    if plugin is None:
        plugin = ReportingPolicyPlugin(config)
        config.pluginmanager.register(plugin, "reporting_policy")
        setattr(config, _REPORTING_POLICY_PLUGIN_ATTR, plugin)

    config.addinivalue_line(
        "markers",
        "pipeline_scenario(name): mark test as covering a key pipeline scenario",
    )

    json_report_file = getattr(config.option, "json_report_file", None)
    if json_report_file:
        try:
            report_path = Path(json_report_file)
        except TypeError:  # pragma: no cover - defensive against unexpected types
            return
        try:
            report_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError:  # pragma: no cover - filesystem issues should not abort tests
            pass
