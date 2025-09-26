from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from library.chembl_client import ChemblClient
from library.config import ApiCfg, MoleculeCatalogCfg
from library.molecule_catalog import fetch_parent_catalog, load_parent_catalog


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


def test_load_parent_catalog_reads_existing_cache(tmp_path: Path, api_cfg: ApiCfg) -> None:
    cache = tmp_path / "catalog.json"
    cache.write_text(json.dumps({"chembl10": "chembl99"}), encoding="utf-8")
    cfg = MoleculeCatalogCfg(cache_path=cache)

    result = load_parent_catalog(
        client=DummyClient([]), api_cfg=api_cfg, catalog_cfg=cfg
    )

    assert result == {"CHEMBL10": "CHEMBL99"}


def test_load_parent_catalog_fetches_and_persists(
    tmp_path: Path, api_cfg: ApiCfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache = tmp_path / "catalog.json"
    cfg = MoleculeCatalogCfg(cache_path=cache)
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
