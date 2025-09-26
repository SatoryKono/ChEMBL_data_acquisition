from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from library.chembl_client import ChemblClient
from library.config import ApiCfg, MoleculeCatalogCfg
from library.molecule_catalog import (
    fetch_parent_catalog,
    fetch_parent_catalog_for,
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
            raise AssertionError("unexpected request")


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
        }
    ]
    client = DummyClient(responses)

    result = fetch_parent_catalog_for(
        [" chembl1 ", "CHEMBL_MISSING"],
        client=client,
        api_cfg=api_cfg,
    )

    assert len(client.calls) == 1
    assert "CHEMBL1%2CCHEMBL_MISSING" in client.calls[0]
    assert result == {"CHEMBL1": "CHEMBL10"}


def test_fetch_parent_catalog_for_chunks_requests(
    api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("library.molecule_catalog._PARENT_LOOKUP_CHUNK_SIZE", 1)
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

    result = fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"],
        client=client,
        api_cfg=api_cfg,
    )

    assert len(client.calls) == 2
    assert result == {"CHEMBL1": "CHEMBL10", "CHEMBL2": "CHEMBL20"}


def test_load_parent_catalog_reads_existing_cache(tmp_path: Path, api_cfg: ApiCfg) -> None:
    cache = tmp_path / "catalog.json"
    cache.write_text(json.dumps({"chembl10": "chembl99"}), encoding="utf-8")
    cfg = MoleculeCatalogCfg(
        cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite"
    )

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
    cfg = MoleculeCatalogCfg(
        cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite"
    )

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == {"CHEMBL1": "CHEMBL42"}


def test_load_parent_catalog_invalid_csv_columns(tmp_path: Path, api_cfg: ApiCfg) -> None:
    cache = tmp_path / "catalog.csv"
    cache.write_text(
        "molecule_chembl_id,parant_molecule_id\nchembl1,chembl42\n",
        encoding="utf-8",
    )
    cfg = MoleculeCatalogCfg(
        cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite"
    )

    with pytest.raises(ValueError):
        load_parent_catalog(client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg)


def test_load_parent_catalog_fetches_and_persists(
    tmp_path: Path, api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache = tmp_path / "catalog.json"
    cfg = MoleculeCatalogCfg(
        cache_path=cache, sqlite_path=tmp_path / "catalog.sqlite"
    )
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
