"""Concurrency regression tests for :mod:`library.molecule_catalog`."""

from __future__ import annotations

import threading
from collections.abc import Mapping

from library.config import ApiCfg, MoleculeCatalogCfg
from library import molecule_catalog


def test_load_parent_catalog_serialises_concurrent_calls(tmp_path, monkeypatch) -> None:
    """Only one thread should fetch the catalogue while others wait for the cache."""

    catalog_data = {"CHEMBL1": "CHEMBL2"}
    fetch_started = threading.Event()
    release_fetch = threading.Event()
    fetch_calls = 0

    def fake_fetch(*args, **kwargs):
        nonlocal fetch_calls
        fetch_calls += 1
        fetch_started.set()
        if not release_fetch.wait(timeout=1):
            raise RuntimeError("fetch not released in time")
        return dict(catalog_data)

    monkeypatch.setattr(molecule_catalog, "fetch_parent_catalog", fake_fetch)

    cfg = MoleculeCatalogCfg(
        cache_path=tmp_path / "catalog.json",
        sqlite_path=tmp_path / "catalog.sqlite",
    )
    api_cfg = ApiCfg(user_agent="chembl-da-tests (mailto:test@example.org)")
    client = object()

    results: dict[str, Mapping[str, str]] = {}
    second_started = threading.Event()

    def first_call() -> None:
        results["first"] = molecule_catalog.load_parent_catalog(
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=cfg,
            timeout=None,
            force_refresh=True,
        )

    def second_call() -> None:
        second_started.set()
        results["second"] = molecule_catalog.load_parent_catalog(
            client=client,
            api_cfg=api_cfg,
            catalog_cfg=cfg,
        )

    first_thread = threading.Thread(target=first_call, name="catalog-first")
    second_thread = threading.Thread(target=second_call, name="catalog-second")

    try:
        first_thread.start()
        assert fetch_started.wait(timeout=1)

        second_thread.start()
        assert second_started.wait(timeout=1)
        assert fetch_calls == 1
        assert second_thread.is_alive()
    finally:
        release_fetch.set()
        first_thread.join(timeout=2)
        second_thread.join(timeout=2)

    assert fetch_calls == 1
    assert results["first"] == catalog_data
    assert results["second"] == catalog_data


def test_helper_skips_catalog_refresh_when_lock_held(monkeypatch) -> None:
    """Helper should not trigger a refresh when the catalogue lock is active."""

    monkeypatch.setattr(
        molecule_catalog,
        "_fetch_parent_catalog_chunk",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(molecule_catalog, "query_parent_catalog", lambda ids, cfg: {})
    monkeypatch.setattr(
        molecule_catalog, "fetch_parent_for_id", lambda *args, **kwargs: None
    )
    monkeypatch.setattr(molecule_catalog, "sleep", lambda delay: None)

    load_calls: list[dict[str, object]] = []

    def fake_load(**kwargs):
        load_calls.append(kwargs)
        return {}

    monkeypatch.setattr(molecule_catalog, "load_parent_catalog", fake_load)

    api_cfg = ApiCfg(user_agent="chembl-da-tests (mailto:test@example.org)")
    cfg = MoleculeCatalogCfg()

    with molecule_catalog._PARENT_CATALOG_LOCK:
        result = molecule_catalog._fetch_parent_catalog_via_helper(
            ["CHEMBL1", "CHEMBL2"],
            client=object(),
            api_cfg=api_cfg,
            timeout=0.1,
            existing={},
            catalog_cfg=cfg,
        )

    assert result == {}
    assert load_calls == []


def test_fetch_parent_catalog_for_preserves_order(monkeypatch) -> None:
    """Fallback single lookups should keep deterministic ordering."""

    api_cfg = ApiCfg(user_agent="chembl-da-tests (mailto:test@example.org)")
    cfg = MoleculeCatalogCfg(filters={})

    monkeypatch.setattr(
        molecule_catalog,
        "_fetch_parent_catalog_chunk",
        lambda *args, **kwargs: {},
    )
    monkeypatch.setattr(molecule_catalog, "query_parent_catalog", lambda *_, **__: {})
    monkeypatch.setattr(molecule_catalog, "load_parent_catalog", lambda **_: {})

    completion_order: list[str] = []
    release_first = threading.Event()

    def fake_fetch_parent_for_id(chembl_id: str, **_: object) -> tuple[str, str]:
        if chembl_id == "CHEMBL1":
            if not release_first.wait(timeout=1):
                raise RuntimeError("second lookup did not complete in time")
            parent = "CHEMBL10"
        else:
            release_first.set()
            parent = "CHEMBL20"
        completion_order.append(chembl_id)
        return chembl_id, parent

    monkeypatch.setattr(
        molecule_catalog, "fetch_parent_for_id", fake_fetch_parent_for_id
    )

    class NoRequestClient:
        def request_json(self, *args, **kwargs):  # pragma: no cover - defensive
            raise AssertionError("unexpected network call")

    client = NoRequestClient()

    result = molecule_catalog.fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"],
        client=client,
        api_cfg=api_cfg,
        catalog_cfg=cfg,
    )

    assert completion_order == ["CHEMBL2", "CHEMBL1"]
    assert list(result.items()) == [
        ("CHEMBL1", "CHEMBL10"),
        ("CHEMBL2", "CHEMBL20"),
    ]
