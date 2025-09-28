from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import requests

from library.clients.chembl_client import ChemblClient
from library.config import ApiCfg, MoleculeCatalogCfg
from library.molecule_catalog import (
    _read_cache,
    fetch_parent_catalog,
    fetch_parent_catalog_for,
    fetch_parent_for_id,
    load_parent_catalog,
    query_parent_catalog,
    update_parent_catalog_cache,
    write_parent_catalog_cache,
)


class DummyClient:
    def __init__(self, responses: list[dict[str, Any]]) -> None:
        self._responses = responses
        self.calls: list[str] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict[str, Any]:
        self.calls.append(url)
        try:
            return self._responses.pop(0)
        except IndexError:  # pragma: no cover - defensive
            raise AssertionError("unexpected request") from None


@pytest.fixture()
def api_cfg() -> ApiCfg:
    return ApiCfg(user_agent="chembl-da-tests (mailto:test@example.org)")


def test_fetch_parent_catalog_normalises_and_paginates(api_cfg: ApiCfg) -> None:
    responses = [
        {
            "molecules": [
                {
                    "molecule_chembl_id": "chembl1",
                    "parent_molecule_chembl_id": "CHEMBL42",
                },
                {
                    "molecule_chembl_id": " chembl2 ",
                    "parent_molecule_chembl_id": " chembl43 ",
                },
            ],
            "page_meta": {"next": "/molecule.json?page=2"},
        },
        {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL3",
                    "parent_molecule_chembl_id": None,
                }
            ],
            "page_meta": {},
        },
    ]
    client = DummyClient(responses)
    cfg = MoleculeCatalogCfg(page_size=2)

    result = fetch_parent_catalog(client=client, api_cfg=api_cfg, catalog_cfg=cfg)

    assert client.calls  # ensure requests were issued
    assert result == {"CHEMBL1": "CHEMBL42", "CHEMBL2": "CHEMBL43"}


def test_fetch_parent_catalog_for_returns_only_requested(api_cfg: ApiCfg) -> None:
    responses = [
        {
            "molecules": [
                {
                    "molecule_chembl_id": "chembl1",
                    "parent_molecule_chembl_id": "chembl10",
                },
                {
                    "molecule_chembl_id": "CHEMBL_EXTRA",
                    "parent_molecule_chembl_id": "CHEMBL99",
                },
            ]
        },
        {"molecule": {"molecule_chembl_id": "CHEMBL_MISSING"}},
    ]
    client = DummyClient(responses)

    result = fetch_parent_catalog_for(
        [" chembl1 ", "CHEMBL_MISSING"],
        client=client,
        api_cfg=api_cfg,
    )

    assert len(client.calls) == 2
    assert "CHEMBL1%2CCHEMBL_MISSING" in client.calls[0]
    assert client.calls[1].endswith(
        "/molecule/CHEMBL_MISSING.json?format=json&fields=molecule_chembl_id%2Cparent_molecule_chembl_id"
    )
    assert result == {"CHEMBL1": "CHEMBL10"}


def test_fetch_parent_catalog_for_chunks_requests(api_cfg: ApiCfg) -> None:
    responses = [
        {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL1",
                    "parent_molecule_chembl_id": "CHEMBL10",
                }
            ]
        },
        {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL2",
                    "parent_molecule_chembl_id": "CHEMBL20",
                }
            ]
        },
    ]
    client = DummyClient(responses)
    cfg = MoleculeCatalogCfg(page_size=1)

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"],
        client=client,
        api_cfg=api_cfg,
        catalog_cfg=cfg,
    )

    assert len(client.calls) == 2
    assert result == {"CHEMBL1": "CHEMBL10", "CHEMBL2": "CHEMBL20"}


def test_fetch_parent_for_id_returns_pair(api_cfg: ApiCfg) -> None:
    responses = [
        {
            "molecule": {
                "molecule_chembl_id": "chembl1",
                "parent_molecule_chembl_id": "chembl10",
            }
        }
    ]
    client = DummyClient(responses)

    pair = fetch_parent_for_id(" chembl1 ", client=client, api_cfg=api_cfg)

    assert pair == ("CHEMBL1", "CHEMBL10")
    assert client.calls == [
        "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL1.json?format=json&fields=molecule_chembl_id%2Cparent_molecule_chembl_id"
    ]


def test_fetch_parent_catalog_for_small_batch_uses_single_helper(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "library.molecule_catalog._PARENT_LOOKUP_FALLBACK_THRESHOLD", 10
    )
    monkeypatch.setattr(
        "library.molecule_catalog.load_parent_catalog", lambda **_: None
    )
    monkeypatch.setattr(
        "library.molecule_catalog.query_parent_catalog", lambda *_, **__: {}
    )
    responses = [
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL1",
                "parent_molecule_chembl_id": "CHEMBL10",
            }
        },
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL2",
                "parent_molecule_chembl_id": "CHEMBL20",
            }
        },
    ]
    client = DummyClient(responses)

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"], client=client, api_cfg=api_cfg
    )

    assert result == {"CHEMBL1": "CHEMBL10", "CHEMBL2": "CHEMBL20"}
    assert all("/molecule/" in call for call in client.calls)


def test_fetch_parent_catalog_for_falls_back_on_bulk_error(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    api_cfg = api_cfg.model_copy(update={"backoff_factor": 0})
    monkeypatch.setattr("library.molecule_catalog._PARENT_LOOKUP_FALLBACK_THRESHOLD", 1)
    monkeypatch.setattr(
        "library.molecule_catalog.load_parent_catalog", lambda **_: None
    )
    monkeypatch.setattr(
        "library.molecule_catalog.query_parent_catalog", lambda *_, **__: {}
    )

    class FallbackClient(DummyClient):
        def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None = None):
            if "molecule.json" in url:
                self.calls.append(url)
                raise requests.RequestException("bulk failure")
            return super().request_json(url, cfg=cfg, timeout=timeout)

    responses = [
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL1",
                "parent_molecule_chembl_id": "CHEMBL10",
            }
        }
    ]
    client = FallbackClient(responses)

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"], client=client, api_cfg=api_cfg
    )

    assert result == {"CHEMBL1": "CHEMBL10"}
    fallback_calls = [call for call in client.calls if "/molecule/" in call]
    assert len(fallback_calls) == api_cfg.retries + 2
    assert any("CHEMBL1" in call for call in fallback_calls)


def test_fetch_parent_catalog_for_partial_fallback_failure(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    api_cfg = api_cfg.model_copy(update={"retries": 2, "backoff_factor": 0})
    monkeypatch.setattr("library.molecule_catalog._PARENT_LOOKUP_FALLBACK_THRESHOLD", 1)

    class PartiallyFailingClient(DummyClient):
        def __init__(self) -> None:
            super().__init__(
                [
                    {
                        "molecule": {
                            "molecule_chembl_id": "CHEMBL1",
                            "parent_molecule_chembl_id": "CHEMBL10",
                        }
                    }
                ]
            )

        def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None = None):
            if "molecule.json" in url:
                self.calls.append(url)
                raise requests.RequestException("bulk failure")
            if "CHEMBL2" in url:
                self.calls.append(url)
                raise requests.RequestException("single failure")
            return super().request_json(url, cfg=cfg, timeout=timeout)

    client = PartiallyFailingClient()

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"], client=client, api_cfg=api_cfg
    )

    assert result == {"CHEMBL1": "CHEMBL10"}
    chembl2_calls = [call for call in client.calls if "/molecule/CHEMBL2" in call]
    assert len(chembl2_calls) == api_cfg.retries + 1


def test_fetch_parent_catalog_for_retries_smaller_batch(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "library.molecule_catalog.load_parent_catalog", lambda **_: None
    )
    monkeypatch.setattr(
        "library.molecule_catalog.query_parent_catalog", lambda *_, **__: {}
    )
    responses = [
        {
            "molecules": [
                {
                    "molecule_chembl_id": "CHEMBL1",
                    "parent_molecule_chembl_id": "CHEMBL10",
                }
            ]
        },
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL2",
                "parent_molecule_chembl_id": "CHEMBL20",
            }
        },
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL3",
                "parent_molecule_chembl_id": "CHEMBL30",
            }
        },
    ]
    client = DummyClient(responses)
    cfg = MoleculeCatalogCfg(page_size=4)

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        client=client,
        api_cfg=api_cfg,
        catalog_cfg=cfg,
    )

    assert result == {
        "CHEMBL1": "CHEMBL10",
        "CHEMBL2": "CHEMBL20",
        "CHEMBL3": "CHEMBL30",
    }
    assert len(client.calls) == 3
    assert "/molecule.json" in client.calls[0]
    assert client.calls[1].endswith(
        "/molecule/CHEMBL2.json?format=json&fields=molecule_chembl_id%2Cparent_molecule_chembl_id"
    )
    assert client.calls[2].endswith(
        "/molecule/CHEMBL3.json?format=json&fields=molecule_chembl_id%2Cparent_molecule_chembl_id"
    )


def test_fetch_parent_catalog_for_respects_single_limit(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "library.molecule_catalog.load_parent_catalog", lambda **_: None
    )
    monkeypatch.setattr(
        "library.molecule_catalog.query_parent_catalog", lambda *_, **__: {}
    )
    responses = [
        {"molecules": []},
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL1",
                "parent_molecule_chembl_id": "CHEMBL10",
            }
        },
        {
            "molecule": {
                "molecule_chembl_id": "CHEMBL2",
                "parent_molecule_chembl_id": "CHEMBL20",
            }
        },
    ]
    client = DummyClient(responses)
    cfg = MoleculeCatalogCfg(page_size=2, fallback_single_limit=1)

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"],
        client=client,
        api_cfg=api_cfg,
        catalog_cfg=cfg,
    )

    assert result == {"CHEMBL1": "CHEMBL10"}
    single_calls = [call for call in client.calls if "/molecule/CHEMBL" in call]
    assert len(single_calls) == 1
    assert "CHEMBL2" not in single_calls[0]


def test_load_parent_catalog_reads_existing_cache(
    tmp_path: Path, api_cfg: ApiCfg
) -> None:
    cache = tmp_path / "catalog.json"
    cache.write_text(json.dumps({"chembl10": "chembl99"}), encoding="utf-8")
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite")

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == {"CHEMBL10": "CHEMBL99"}


def test_load_parent_catalog_reads_csv_cache(tmp_path: Path, api_cfg: ApiCfg) -> None:
    cache = tmp_path / "catalog.csv"
    cache.write_text(
        "molecule_chembl_id,parent_molecule_chembl_id\nchembl1,chembl42\n",
        encoding="utf-8",
    )
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite")

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == {"CHEMBL1": "CHEMBL42"}


def test_load_parent_catalog_invalid_csv_columns(
    tmp_path: Path, api_cfg: ApiCfg
) -> None:
    cache = tmp_path / "catalog.csv"
    cache.write_text(
        "molecule_chembl_id,parant_molecule_id\nchembl1,chembl42\n",
        encoding="utf-8",
    )
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite")

    with pytest.raises(ValueError):
        load_parent_catalog(client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg)


def test_load_parent_catalog_fetches_and_persists(
    tmp_path: Path, api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache = tmp_path / "catalog.json"
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite")
    data = {"CHEMBL50": "CHEMBL60"}

    def fake_fetch(
        *,
        client: ChemblClient,
        api_cfg: ApiCfg,
        catalog_cfg: MoleculeCatalogCfg,
        timeout: float | None = None,
    ) -> dict[str, str]:
        return data

    monkeypatch.setattr(
        "library.molecule_catalog.fetch_parent_catalog",
        fake_fetch,
    )

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == data
    assert json.loads(cache.read_text(encoding="utf-8")) == data
    assert cfg.sqlite_path.is_file()


def test_query_parent_catalog_reads_sqlite(tmp_path: Path) -> None:
    cache = tmp_path / "catalog.json"
    sqlite_path = tmp_path / "catalog.sqlite"
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=sqlite_path)
    write_parent_catalog_cache({"CHEMBL1": "CHEMBL10", "CHEMBL2": "CHEMBL20"}, cfg)

    result = query_parent_catalog(["chembl1", "chembl3"], cfg)

    assert result == {"CHEMBL1": "CHEMBL10"}
    assert sqlite_path.is_file()


def test_update_parent_catalog_cache_appends(tmp_path: Path, api_cfg: ApiCfg) -> None:
    cache = tmp_path / "catalog.json"
    sqlite_path = tmp_path / "catalog.sqlite"
    cfg = MoleculeCatalogCfg(cache_path=cache, sqlite_path=sqlite_path)
    write_parent_catalog_cache({"CHEMBL1": "CHEMBL10"}, cfg)

    update_parent_catalog_cache({"CHEMBL2": "CHEMBL20"}, cfg)

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == {"CHEMBL1": "CHEMBL10", "CHEMBL2": "CHEMBL20"}
    assert json.loads(cache.read_text(encoding="utf-8")) == result


def test_read_cache_invalid_json_returns_empty(tmp_path: Path) -> None:
    """Gracefully ignore cache files containing malformed JSON."""

    cache_path = tmp_path / "cache.json"
    cache_path.write_text("[1, 2, 3]", encoding="utf-8")
    cfg = MoleculeCatalogCfg(cache_path=cache_path)

    assert _read_cache(cache_path, cfg) == {}


def test_read_cache_missing_csv_columns_raises(tmp_path: Path) -> None:
    """CSV caches missing required headers should raise an explicit error."""

    cache_path = tmp_path / "cache.csv"
    cache_path.write_text("child\nCHEMBL1\n", encoding="utf-8")
    cfg = MoleculeCatalogCfg(
        cache_path=cache_path, child_field="child", parent_field="parent"
    )

    with pytest.raises(ValueError, match="missing columns"):
        _read_cache(cache_path, cfg)
