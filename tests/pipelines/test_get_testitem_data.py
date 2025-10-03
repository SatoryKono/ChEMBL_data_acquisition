from __future__ import annotations

import argparse
import json

import threading
import time

from collections.abc import (
    Callable,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from typing import Any
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import requests

from library.integration import chembl_library as cl
from library import io
from library.integration import molecule_catalog
from library.integration import pubchem_library as pl
import library.testitem_pipeline as pipeline
import library.testitem_pipeline.cli as pipeline_cli
import library.testitem_pipeline.catalog as modern_catalog_module
import library.pipelines.testitem as legacy_pipeline

from library.config import ApiCfg, Config, IoCfg, MoleculeCatalogCfg

from library.schemas import TestitemsSchema
from scripts import get_testitem_data as gtd


def test_ensure_no_parant_column_smoke() -> None:
    """The helper re-exports the pipeline validation logic."""

    valid_df = pd.DataFrame({"parent_molecule_id": ["CHEMBL1"]})
    gtd.ensure_no_parant_column(valid_df)

    invalid_df = pd.DataFrame({pipeline._TYPO_PARENT_COLUMN: ["CHEMBL1"]})
    with pytest.raises(ValueError):
        gtd.ensure_no_parant_column(invalid_df)


def prepare_parent_lookup_data(
    df: pd.DataFrame, catalog_cfg
) -> pipeline.ParentLookupPreparedData:
    child_column = catalog_cfg.child_field
    parent_column = catalog_cfg.parent_field

    if child_column in df.columns:
        normalised_child = (
            df[child_column].astype("string").fillna("").str.strip().str.upper()
        )
    else:
        normalised_child = pd.Series("", index=df.index, dtype="string")

    if parent_column in df.columns:
        existing_parent = (
            df[parent_column].astype("string").fillna("").str.strip().str.upper()
        )
    else:
        existing_parent = pd.Series("", index=df.index, dtype="string")

    need_lookup_mask = (normalised_child != "") & (existing_parent == "")
    need_lookup = set(normalised_child[need_lookup_mask])

    return pipeline.ParentLookupPreparedData(
        child_ids=normalised_child,
        existing_parent_ids=existing_parent,
        need_lookup=need_lookup,
    )


def test_read_input_ids_limit_and_sample(tmp_path: Path, cfg: Config) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text(
        "molecule_chembl_id\nCHEMBL1\nCHEMBL2\n",
        encoding=cfg.io.csv_encoding,
    )

    status, result = gtd.read_input_ids(
        input_csv,
        column=cfg.testitem.column,
        io_cfg=cfg.io,
        limit=1,
    )

    assert status == 0
    assert result is not None
    assert list(result.ids_iter) == ["CHEMBL1"]
    assert result.sample_ids == ("CHEMBL1",)


def test_read_input_ids_offset(tmp_path: Path, cfg: Config) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text(
        "molecule_chembl_id\nCHEMBL1\nCHEMBL2\nCHEMBL3\n",
        encoding=cfg.io.csv_encoding,
    )

    status, result = gtd.read_input_ids(
        input_csv,
        column=cfg.testitem.column,
        io_cfg=cfg.io,
        limit=None,
        offset=2,
    )

    assert status == 0
    assert result is not None
    assert list(result.ids_iter) == ["CHEMBL3"]
    assert result.sample_ids == ("CHEMBL3",)


def test_read_input_ids_missing_file(tmp_path: Path, cfg: Config) -> None:
    status, result = gtd.read_input_ids(
        tmp_path / "missing.csv",
        column=cfg.testitem.column,
        io_cfg=cfg.io,
        limit=None,
    )

    assert status == 1
    assert result is None


def test_fetch_testitems_failure(monkeypatch: pytest.MonkeyPatch, cfg: Config) -> None:
    def fail_fetch(*args: object, **kwargs: object) -> pd.DataFrame:
        raise requests.RequestException("boom")

    monkeypatch.setattr(cl, "get_testitem", fail_fetch)

    status, chunks, requested_ids = gtd.fetch_testitems(
        iter(["CHEMBL1"]),
        api_cfg=cfg.api,
        batch_size=1,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        sample_ids=("CHEMBL1",),
        fields=cfg.testitem.fields,
        page_limit=cfg.testitem.request_limit,
        retry_cfg=cfg.testitem.batch_retry,
    )

    assert status == 1
    assert chunks is None
    assert requested_ids == ()


def test_fetch_testitems_retries_reduced_batch(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    attempts: list[int] = []

    def flaky_fetch(
        ids: Iterable[str],
        *,
        cfg: ApiCfg,
        client: object,
        chunk_size: int,
        timeout: float | None,
        fields: Sequence[str] | None = None,
        page_limit: int = 0,
    ) -> pd.DataFrame:
        _ = cfg, client, timeout, fields, page_limit
        attempts.append(chunk_size)
        if chunk_size >= 4:
            raise requests.RequestException("boom")
        values = list(ids)
        return pd.DataFrame({"molecule_chembl_id": values})

    monkeypatch.setattr(cl, "get_testitem", flaky_fetch)

    warnings: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, *args: object, **kwargs: object) -> None:
        warnings.append((event, dict(kwargs)))

    monkeypatch.setattr(gtd.logger, "warning", capture_warning)

    original_retry = cfg.testitem.batch_retry
    cfg.testitem.batch_retry = cfg.testitem.batch_retry.model_copy(
        update={"enable": True, "min_size": 1, "shrink_factor": 0.5}
    )

    status, chunks, requested_ids = gtd.fetch_testitems(
        iter(["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4"]),
        api_cfg=cfg.api,
        batch_size=4,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        sample_ids=("CHEMBL1",),
        fields=cfg.testitem.fields,
        page_limit=cfg.testitem.request_limit,
        retry_cfg=cfg.testitem.batch_retry,
    )

    cfg.testitem.batch_retry = original_retry

    assert status == 0
    assert chunks is not None
    frames = list(chunks)
    assert len(frames) == 1
    assert frames[0]["molecule_chembl_id"].tolist() == [
        "CHEMBL1",
        "CHEMBL2",
        "CHEMBL3",
        "CHEMBL4",
    ]
    assert attempts == [4, 2]
    assert requested_ids == (
        "CHEMBL1",
        "CHEMBL2",
        "CHEMBL3",
        "CHEMBL4",
    )
    retry_events = [event for event, _ in warnings if event == "testitem_fetch_retry_reduced_batch"]
    assert retry_events


def test_ensure_no_parant_column_raises() -> None:
    df = pd.DataFrame([{"parant_molecule_id": "CHEMBL1"}])

    with pytest.raises(
        ValueError,
        match="unexpected column 'parant_molecule_id'; use 'parent_molecule_id' instead",
    ):
        gtd.ensure_no_parant_column(df)


def test_fetch_testitems_passes_fields_and_limit(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    captured: dict[str, object] = {}

    def fake_get_testitem(
        ids: Iterable[str],
        *,
        cfg: ApiCfg,
        client: object,
        chunk_size: int,
        timeout: float | None,
        fields: Sequence[str] | None = None,
        page_limit: int = 0,
    ) -> pd.DataFrame:
        captured["fields"] = fields
        captured["page_limit"] = page_limit
        values = list(ids)
        return pd.DataFrame(
            [{"molecule_chembl_id": value} for value in values],
            columns=["molecule_chembl_id"],
        )

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)

    status, chunks, requested_ids = gtd.fetch_testitems(
        iter(["CHEMBL1", "CHEMBL2"]),
        api_cfg=cfg.api,
        batch_size=2,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        sample_ids=("CHEMBL1",),
        fields=("a", "b"),
        page_limit=500,
        retry_cfg=cfg.testitem.batch_retry,
    )

    assert status == 0
    assert chunks is not None
    frames = list(chunks)
    assert len(frames) == 1
    assert frames[0]["molecule_chembl_id"].tolist() == ["CHEMBL1", "CHEMBL2"]
    assert captured["fields"] == ("a", "b")
    assert captured["page_limit"] == 500
    assert requested_ids == ("CHEMBL1", "CHEMBL2")


def test_fetch_testitems_logs_missing_summary(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    captured: list[tuple[str, dict[str, object]]] = []

    def fake_warning(event: str, *args: object, **kwargs: object) -> None:
        captured.append((event, dict(kwargs)))

    monkeypatch.setattr(gtd.logger, "warning", fake_warning)

    def fake_get_testitem(
        ids: Iterable[str],
        *,
        cfg: ApiCfg,
        client: object,
        chunk_size: int,
        timeout: float | None,
        fields: Sequence[str] | None = None,
        page_limit: int = 0,
    ) -> pd.DataFrame:
        _ = list(ids)
        return pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)

    status, chunks, requested_ids = gtd.fetch_testitems(
        iter(["CHEMBL1", "chembl2", "CHEMBL3"]),
        api_cfg=cfg.api,
        batch_size=2,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        sample_ids=("CHEMBL1",),
        fields=cfg.testitem.fields,
        page_limit=500,
        retry_cfg=cfg.testitem.batch_retry,
    )

    assert status == 0
    assert chunks is not None
    assert requested_ids == ("CHEMBL1", "chembl2", "CHEMBL3")
    _ = list(chunks)
    missing = [
        record for record in captured if record[0] == "chembl_missing_identifiers"
    ]
    assert not missing


def test_log_missing_identifier_summary_limits_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: list[tuple[str, dict[str, object]]] = []

    def fake_warning(event: str, *args: object, **kwargs: object) -> None:
        captured.append((event, dict(kwargs)))

    monkeypatch.setattr(pipeline.logger, "warning", fake_warning)

    pipeline._log_missing_identifier_summary([])
    assert captured == []

    identifiers = [f"CHEMBL{index}" for index in range(20)]
    pipeline._log_missing_identifier_summary(identifiers)

    assert len(captured) == 1
    event, payload = captured[0]
    assert event == "chembl_missing_identifiers"
    assert payload["total"] == len(identifiers)
    assert payload["sample"] == identifiers[:10]
    assert "identifiers" not in payload


def test_fetch_parent_catalog_skips_single_when_parentless(
    cfg: Config,
) -> None:
    captured_urls: list[str] = []

    class DummyClient:
        def request_json(
            self,
            url: str,
            *,
            cfg: ApiCfg,
            timeout: float | None,
        ) -> Mapping[str, object]:
            captured_urls.append(url)
            if "/molecule.json" in url:
                return {"molecules": []}
            raise AssertionError("single molecule fetch should be skipped")

    client = DummyClient()
    cfg.molecule_catalog.filters = {"parent_molecule_chembl_id__isnull": "false"}

    result = molecule_catalog.fetch_parent_catalog_for(
        ["CHEMBL1", "CHEMBL2"],
        client=client,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        catalog_cfg=cfg.molecule_catalog,
    )

    assert result == {}
    assert captured_urls
    assert all("/molecule.json" in url for url in captured_urls)


def test_prepare_pubchem_api_cfg_prefers_pubchem_override(cfg: Config) -> None:
    cfg.pubchem.user_agent = "chembl-da/1.0 (mailto:pubchem@example.org)"

    pubchem_cfg = gtd._prepare_pubchem_api_cfg(cfg, cfg.api)

    assert pubchem_cfg.user_agent == cfg.pubchem.user_agent
    assert pubchem_cfg is not cfg.api


def test_prepare_pubchem_api_cfg_reexport_matches_pipeline() -> None:
    assert gtd._prepare_pubchem_api_cfg is pipeline._prepare_pubchem_api_cfg


def test_prepare_pubchem_api_cfg_requires_custom_user_agent(cfg: Config) -> None:
    placeholder = "chembl-da/0.1 (mailto:contact@example.org)"
    cfg.api.user_agent = placeholder
    cfg.pubchem.user_agent = placeholder

    with pytest.raises(ValueError, match="PubChem configuration requires a user_agent"):
        gtd._prepare_pubchem_api_cfg(cfg, cfg.api)


def test_fetch_parent_catalog_single_entry_parentless_uses_bulk_only(
    cfg: Config,
) -> None:
    captured_urls: list[str] = []

    class DummyClient:
        def request_json(
            self,
            url: str,
            *,
            cfg: ApiCfg,
            timeout: float | None,
        ) -> Mapping[str, object]:
            captured_urls.append(url)
            if "/molecule.json" in url:
                return {"molecules": []}
            raise AssertionError("single molecule fetch should be skipped")

    client = DummyClient()
    cfg.molecule_catalog.filters = {"parent_molecule_chembl_id__isnull": "false"}

    result = molecule_catalog.fetch_parent_catalog_for(
        ["CHEMBL1"],
        client=client,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        catalog_cfg=cfg.molecule_catalog,
    )

    assert result == {}
    assert captured_urls
    assert all("/molecule.json" in url for url in captured_urls)


def test_prepare_parent_enrichment_uses_lookup_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    df = pd.DataFrame({child_field: ["CHEMBL1"], parent_field: [pd.NA]})
    cfg.molecule_catalog.cache_path = tmp_path / "cache.json"
    cfg.molecule_catalog.sqlite_path = tmp_path / "cache.sqlite"
    hierarchy_path = tmp_path / "hierarchy.csv"

    captured_path: Path | None = None

    def fake_lookup(path: Path | None, *, io_cfg: IoCfg, **_: object) -> dict[str, str]:
        nonlocal captured_path
        captured_path = path
        return {"CHEMBL1": "CHEMBL999"}

    for target in (pipeline, pipeline.catalog):
        monkeypatch.setattr(target, "load_molecule_hierarchy_lookup", fake_lookup)
    monkeypatch.setattr(pipeline, "query_parent_catalog", lambda *_, **__: {})
    monkeypatch.setattr(
        pipeline.molecule_catalog, "fetch_parent_catalog_for", lambda *_, **__: {}
    )
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda *_, **__: {})

    status, prep = pipeline.prepare_parent_enrichment(
        df.copy(),
        catalog_cfg=cfg.molecule_catalog,
        io_cfg=cfg.io,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        hierarchy_lookup_path=hierarchy_path,
    )

    assert status == 0
    assert prep is not None
    assert captured_path == hierarchy_path
    assert list(prep.lookup_data.child_ids) == ["CHEMBL1"]


def test_prepare_parent_enrichment_dictionary_sets_parent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    df = pd.DataFrame({child_field: ["CHEMBL1"], parent_field: [pd.NA]})
    cfg.molecule_catalog.cache_path = tmp_path / "cache.json"
    cfg.molecule_catalog.cache_path.write_text("{}")
    cfg.molecule_catalog.sqlite_path = tmp_path / "cache.sqlite"
    hierarchy_path = tmp_path / "hierarchy.csv"
    hierarchy_path.write_text(
        "molecule_chembl_id,parent_molecule_chembl_id\nCHEMBL1,CHEMBL999\n",
        encoding=cfg.io.csv_encoding,
    )
    cfg.molecule_catalog.hierarchy_lookup_path = hierarchy_path
    cfg.molecule_catalog.hierarchy_lookup_encoding = cfg.io.csv_encoding
    cfg.molecule_catalog.hierarchy_lookup_delimiter = cfg.io.csv_sep
    for target in (pipeline, pipeline.catalog):
        monkeypatch.setattr(target, "query_parent_catalog", lambda *_, **__: pytest.fail("no fallback"))
        monkeypatch.setattr(target, "load_parent_catalog", lambda *_, **__: {})
        monkeypatch.setattr(target, "update_parent_catalog_cache", lambda *_, **__: None)
        monkeypatch.setattr(target, "write_parent_catalog_cache", lambda *_, **__: None)

    status, prep = pipeline.prepare_parent_enrichment(
        df.copy(),
        catalog_cfg=cfg.molecule_catalog,
        io_cfg=cfg.io,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        hierarchy_lookup_path=None,
    )

    assert status == 0
    assert prep is not None
    assert prep.lookup_data.need_lookup == set()
    assert prep.parent_stats.missing == 0
    assert prep.parent_stats.attached == 1
    assert prep.df[parent_field].tolist() == ["CHEMBL999"]


def test_prepare_parent_enrichment_dictionary_null_parent_skips_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    df = pd.DataFrame({child_field: ["CHEMBL1"], parent_field: [pd.NA]})
    cfg.molecule_catalog.cache_path = tmp_path / "cache.json"
    cfg.molecule_catalog.cache_path.write_text("{}")
    cfg.molecule_catalog.sqlite_path = tmp_path / "cache.sqlite"
    hierarchy_path = tmp_path / "hierarchy.csv"
    hierarchy_path.write_text(
        "molecule_chembl_id,parent_molecule_chembl_id\nCHEMBL1,NULL\n",
        encoding=cfg.io.csv_encoding,
    )
    cfg.molecule_catalog.hierarchy_lookup_path = hierarchy_path
    cfg.molecule_catalog.hierarchy_lookup_encoding = cfg.io.csv_encoding
    cfg.molecule_catalog.hierarchy_lookup_delimiter = cfg.io.csv_sep
    for target in (pipeline, pipeline.catalog):
        monkeypatch.setattr(target, "query_parent_catalog", lambda *_, **__: pytest.fail("no fallback"))
        monkeypatch.setattr(target, "load_parent_catalog", lambda *_, **__: {})
        monkeypatch.setattr(target, "update_parent_catalog_cache", lambda *_, **__: None)
        monkeypatch.setattr(target, "write_parent_catalog_cache", lambda *_, **__: None)

    status, prep = pipeline.prepare_parent_enrichment(
        df.copy(),
        catalog_cfg=cfg.molecule_catalog,
        io_cfg=cfg.io,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        hierarchy_lookup_path=None,
    )

    assert status == 0
    assert prep is not None
    assert prep.lookup_data.need_lookup == set()
    assert prep.parent_stats.missing == 0
    assert prep.parent_stats.attached == 1
    assert pd.isna(prep.df[parent_field].iloc[0])


def test_prepare_parent_enrichment_falls_back_when_missing_from_dictionary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    df = pd.DataFrame({child_field: ["CHEMBL1"], parent_field: [pd.NA]})
    cfg.molecule_catalog.cache_path = tmp_path / "cache.json"
    cfg.molecule_catalog.cache_path.write_text("{}")
    cfg.molecule_catalog.sqlite_path = tmp_path / "cache.sqlite"
    hierarchy_path = tmp_path / "hierarchy.csv"
    hierarchy_path.write_text(
        "molecule_chembl_id,parent_molecule_chembl_id\nCHEMBL2,CHEMBL3\n",
        encoding=cfg.io.csv_encoding,
    )
    cfg.molecule_catalog.hierarchy_lookup_path = hierarchy_path
    cfg.molecule_catalog.hierarchy_lookup_encoding = cfg.io.csv_encoding
    cfg.molecule_catalog.hierarchy_lookup_delimiter = cfg.io.csv_sep

    captured: dict[str, set[str]] = {"need": set()}

    def fake_query(values: set[str], *, catalog_cfg: MoleculeCatalogCfg) -> dict[str, str]:
        captured["need"] = set(values)
        return {}

    for target in (pipeline, pipeline.catalog):
        monkeypatch.setattr(target, "query_parent_catalog", fake_query)
        monkeypatch.setattr(target, "load_parent_catalog", lambda *_, **__: {})
        monkeypatch.setattr(target, "update_parent_catalog_cache", lambda *_, **__: None)
        monkeypatch.setattr(target, "write_parent_catalog_cache", lambda *_, **__: None)

    status, prep = pipeline.prepare_parent_enrichment(
        df.copy(),
        catalog_cfg=cfg.molecule_catalog,
        io_cfg=cfg.io,
        api_cfg=cfg.api,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        hierarchy_lookup_path=None,
    )

    assert status == 0
    assert prep is not None
    assert prep.lookup_data.need_lookup == {"CHEMBL1"}
    assert captured["need"] == {"CHEMBL1"}
    assert prep.parent_stats.missing == 1
    assert prep.parent_stats.attached == 0


def test_run_parent_enrichment_failure(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({cfg.molecule_catalog.parent_field: [pd.NA]})
    lookup = pipeline.ParentLookupPreparedData(
        child_ids=pd.Series(["CHEMBL1"], dtype="string"),
        existing_parent_ids=pd.Series([""], dtype="string"),
        need_lookup=set(),
    )
    prep = pipeline.ParentEnrichmentPreparation(
        df=df,
        lookup_data=lookup,
        parent_catalog=None,
        parent_catalog_source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        parent_stats=pipeline.ParentLookupStats(
            source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
            missing=0,
            unique=0,
            attached=0,
            uncovered=0,
        ),
    )

    def fail_attach(*args: object, **kwargs: object) -> tuple[pd.DataFrame, object]:
        raise ValueError("attach failed")

    monkeypatch.setattr(pipeline, "attach_parent_molecule_ids", fail_attach)

    status, result = pipeline.run_parent_enrichment(
        prep,
        client=SimpleNamespace(),
        api_cfg=cfg.api,
        catalog_cfg=cfg.molecule_catalog,
        timeout=cfg.testitem.timeout,
    )

    assert status == 1
    assert result is None


def test_augment_pubchem_initialises_caches(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    cfg.pubchem.enable = True
    cfg.pubchem.cid_cache_path = tmp_path / "pubchem_cache.json"
    cfg.pubchem.cache_ttl_hours = 24

    captured: dict[str, object] = {}

    def fake_load(
        path: Path | None, ttl_hours: float | None = None
    ) -> dict[str, str | None]:
        captured["load_args"] = (path, ttl_hours)
        return {}

    def fake_add(
        frame: pd.DataFrame,
        pubchem_cfg: pl.PubChemCfg,
        **kwargs: object,
    ) -> pd.DataFrame:
        captured["add_kwargs"] = kwargs
        return frame.assign(pubchem_cid="1")

    monkeypatch.setattr(pipeline, "_load_pubchem_cid_cache", fake_load)
    monkeypatch.setattr(pipeline, "add_pubchem_data", fake_add)

    pubchem_api_cfg = pipeline._prepare_pubchem_api_cfg(cfg, cfg.api)

    result = pipeline.augment_pubchem(
        df,
        pubchem_cfg=cfg.pubchem,
        api_cfg=pubchem_api_cfg,
        retry_cfg=cfg.retry,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        fields=cfg.testitem.fields,
        request_limit=cfg.testitem.request_limit,
    )

    assert captured["load_args"] == (
        cfg.pubchem.cid_cache_path,
        cfg.pubchem.cache_ttl_hours,
    )
    assert "cid_cache" in captured["add_kwargs"]
    assert captured["add_kwargs"].get("testitem_fields") == cfg.testitem.fields
    assert captured["add_kwargs"].get("request_limit") == cfg.testitem.request_limit
    assert "pubchem_cid" in result.columns


def test_augment_pubchem_initialises_session_once(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    cfg.pubchem.enable = True

    monkeypatch.setattr(pipeline, "_PUBCHEM_SESSION_SIGNATURE", None, raising=False)

    captured: dict[str, object] = {"init_calls": 0}

    pubchem_api_cfg = pipeline._prepare_pubchem_api_cfg(cfg, cfg.api)

    def fake_init_session(api: object, retry: object) -> None:
        captured["init_calls"] = captured["init_calls"] + 1
        captured["init_args"] = (api, retry)

    def fake_add(
        frame: pd.DataFrame,
        pubchem_cfg: pl.PubChemCfg,
        **kwargs: object,
    ) -> pd.DataFrame:
        return frame.assign(pubchem_cid="1")

    monkeypatch.setattr(pipeline.pl, "init_session", fake_init_session)
    monkeypatch.setattr(pipeline, "add_pubchem_data", fake_add)

    first = pipeline.augment_pubchem(
        df,
        pubchem_cfg=cfg.pubchem,
        api_cfg=pubchem_api_cfg,
        retry_cfg=cfg.retry,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        fields=cfg.testitem.fields,
        request_limit=cfg.testitem.request_limit,
    )

    second = pipeline.augment_pubchem(
        df,
        pubchem_cfg=cfg.pubchem,
        api_cfg=pubchem_api_cfg,
        retry_cfg=cfg.retry,
        timeout=cfg.testitem.timeout,
        client=SimpleNamespace(),
        fields=cfg.testitem.fields,
        request_limit=cfg.testitem.request_limit,
    )

    assert captured["init_calls"] == 1
    assert captured["init_args"] == (pubchem_api_cfg, cfg.retry)
    assert "pubchem_cid" in first.columns
    assert "pubchem_cid" in second.columns


def test_apply_testitem_enrichment_disable(cfg: Config) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    status, result = pipeline.apply_testitem_enrichment(
        df,
        enrichment_cfg=SimpleNamespace(enable=False),
        io_cfg=cfg.io,
    )

    assert status == 0
    assert result is df


def test_apply_testitem_enrichment_failure(
    monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    def fail_enrich(*args: object, **kwargs: object) -> pd.DataFrame:
        raise ValueError("enrich failed")

    monkeypatch.setattr(pipeline.testitem_enrichment, "enrich", fail_enrich)

    status, result = pipeline.apply_testitem_enrichment(
        df,
        enrichment_cfg=SimpleNamespace(enable=True),
        io_cfg=cfg.io,
    )

    assert status == 1
    assert result is None


def test_finalize_output_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    parent_stats = pipeline.ParentLookupStats(
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        missing=0,
        unique=1,
        attached=1,
        uncovered=0,
        failed_ids=("CHEMBL9",),
    )

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> SimpleNamespace:
        return SimpleNamespace(data=frame, failure_cases=pd.DataFrame())

    monkeypatch.setattr(pipeline, "validate_testitems", fake_validate, raising=False)
    captured_stats: dict[str, object] = {}

    def capture_meta(**kwargs: object) -> None:
        stats = kwargs.get("stats")
        if isinstance(stats, dict):
            captured_stats.update(stats)

    monkeypatch.setattr(pipeline, "write_meta_yaml", capture_meta)
    monkeypatch.setattr(pipeline_cli, "write_meta_yaml", capture_meta)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(pipeline_cli, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(
        pipeline, "analyze_table_quality", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(
        pipeline_cli, "analyze_table_quality", lambda *_args, **_kwargs: None
    )

    written_chunks: list[pd.DataFrame] = []

    def capture_chunks(
        chunks: Iterable[pd.DataFrame],
        path: Path,
        **kwargs: object,
    ) -> Path:
        materialized = list(chunks)
        written_chunks.extend(materialized)
        return path

    monkeypatch.setattr(
        pipeline,
        "write_csv_chunks_deterministic",
        capture_chunks,
        raising=False,
    )
    monkeypatch.setattr(
        pipeline_cli,
        "write_csv_chunks_deterministic",
        capture_chunks,
    )

    exit_code = pipeline.finalize_output(
        [df],
        cfg=cfg,
        output=tmp_path / "out.csv",
        parent_stats_supplier=lambda: parent_stats,
        input_csv=tmp_path / "in.csv",
    )

    assert exit_code == 0
    assert captured_stats.get("parent_lookup_failed_ids") == ["CHEMBL9"]
    assert captured_stats.get("parent_lookup_failed_count") == 1

    written_columns: set[str] = set()
    for chunk in written_chunks:
        written_columns.update(chunk.columns)

    assert written_columns == set(TestitemsSchema.columns)


def test_finalize_output_missing_required_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"extra": ["value"]})
    parent_stats = pipeline.ParentLookupStats(
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )

    monkeypatch.setattr(
        pipeline,
        "write_csv_chunks_deterministic",
        lambda *args, **kwargs: pytest.fail(
            "should not write output when required columns are missing"
        ),
        raising=False,
    )
    monkeypatch.setattr(
        pipeline,
        "analyze_table_quality",
        lambda *args, **kwargs: pytest.fail(
            "should not analyze quality when required columns are missing"
        ),
    )

    exit_code = pipeline.finalize_output(
        [df],
        cfg=cfg,
        output=tmp_path / "out.csv",
        parent_stats_supplier=lambda: parent_stats,
        input_csv=tmp_path / "in.csv",
    )

    assert exit_code == 1
    assert not (tmp_path / "out.csv").exists()


def test_finalize_output_validation_failure_skips_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    parent_stats = pipeline.ParentLookupStats(
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        missing=0,
        unique=1,
        attached=1,
        uncovered=0,
    )

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> SimpleNamespace:
        failures = pd.DataFrame([{"column": "value", "failure": "bad"}])
        return SimpleNamespace(data=frame, failure_cases=failures)

    monkeypatch.setattr(pipeline_cli, "validate_testitems", fake_validate)
    monkeypatch.setattr(
        pipeline_cli,
        "write_csv_chunks_deterministic",
        lambda *args, **kwargs: pytest.fail("should not write output when validation fails"),
    )
    monkeypatch.setattr(
        pipeline_cli,
        "write_meta_yaml",
        lambda *args, **kwargs: pytest.fail("should not write metadata when validation fails"),
    )
    monkeypatch.setattr(
        pipeline_cli,
        "file_sha256",
        lambda *args, **kwargs: pytest.fail("should not hash output when validation fails"),
    )
    monkeypatch.setattr(
        pipeline_cli,
        "analyze_table_quality",
        lambda *args, **kwargs: pytest.fail("should not analyze quality when validation fails"),
    )

    output_path = tmp_path / "out.csv"
    exit_code = pipeline.finalize_output(
        [df],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=lambda: parent_stats,
        input_csv=tmp_path / "in.csv",
    )

    assert exit_code == 1
    assert not output_path.exists()
    assert not (output_path.with_name(output_path.name + ".meta.yaml")).exists()
    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
    assert failure_path.exists()
    assert (failure_path.with_name(failure_path.name + ".meta.yaml")).exists()


def test_finalize_output_streams_sorted_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": [
                "CHEMBL3",
                "CHEMBL1",
                "CHEMBL2",
                "CHEMBL5",
            ],
            "value": [3, 1, 2, 4],
        }
    )
    parent_stats = pipeline.ParentLookupStats(
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        missing=0,
        unique=4,
        attached=4,
        uncovered=0,
    )

    def fake_validate(frame: pd.DataFrame, *, return_result: bool) -> SimpleNamespace:
        return SimpleNamespace(data=frame, failure_cases=pd.DataFrame())

    recorded_targets: list[str] = []
    original_to_csv = pd.DataFrame.to_csv

    def spy_to_csv(self, *args, **kwargs):  # type: ignore[override]
        target = args[0] if args else kwargs.get("path_or_buf")
        recorded_targets.append(str(target))
        return original_to_csv(self, *args, **kwargs)

    cfg.io.csv_chunksize = 2

    monkeypatch.setattr(pipeline, "validate_testitems", fake_validate)
    monkeypatch.setattr(pd.DataFrame, "to_csv", spy_to_csv)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "hash")
    monkeypatch.setattr(
        pipeline, "analyze_table_quality", lambda df, table_name, **_: None
    )

    exit_code = pipeline.finalize_output(
        [df],
        cfg=cfg,
        output=tmp_path / "out.csv",
        parent_stats_supplier=lambda: parent_stats,
        input_csv=tmp_path / "in.csv",
    )

    chunk_writes = [target for target in recorded_targets if "chunk_" in target]
    assert len(chunk_writes) >= 2
    assert exit_code == 0


def test_load_molecule_hierarchy_lookup_missing(tmp_path: Path, cfg: Config) -> None:
    path = tmp_path / "missing.csv"

    result = pipeline.load_molecule_hierarchy_lookup(path, io_cfg=cfg.io)

    assert result == {}


def test_load_molecule_hierarchy_lookup_filters_empty_rows(
    tmp_path: Path, cfg: Config
) -> None:
    path = tmp_path / "hierarchy.csv"
    path.write_text(
        "\n".join(
            [
                "molecule_chembl_id,parent_molecule_chembl_id",
                "CHEMBL1,CHEMBL2",
                "CHEMBL3,",
                "CHEMBL4,NULL",
                ",CHEMBL4",
                " chembl5 , chembl6 ",
                "CHEMBL1,CHEMBL7",
            ]
        ),
        encoding=cfg.io.csv_encoding,
    )

    result = pipeline.load_molecule_hierarchy_lookup(path, io_cfg=cfg.io)

    assert result == {
        "CHEMBL1": "CHEMBL2",
        "CHEMBL3": None,
        "CHEMBL4": None,
        "CHEMBL5": "CHEMBL6",
    }


def test_load_molecule_hierarchy_lookup_missing_columns(
    tmp_path: Path, cfg: Config, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "hierarchy_invalid.csv"
    path.write_text(
        "molecule_chembl_id,parent_id\nCHEMBL1,CHEMBL2\n",
        encoding=cfg.io.csv_encoding,
    )

    captured: list[tuple[str, dict[str, object]]] = []

    def fake_warning(event: str, *args: object, **kwargs: object) -> None:
        captured.append((event, dict(kwargs)))

    monkeypatch.setattr(pipeline.logger, "warning", fake_warning)

    result = pipeline.load_molecule_hierarchy_lookup(path, io_cfg=cfg.io)

    assert result == {"CHEMBL1": None}
    assert captured
    event, payload = captured[-1]
    assert event == "molecule_hierarchy_missing_parent_column"
    assert payload.get("column") == "parent_molecule_chembl_id"


def test_prepare_pubchem_caches_primes_local_parent_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"],
            "parent_molecule_chembl_id": [pd.NA, pd.NA, "CHEMBL1"],
            "canonical_smiles": ["C", "N", "CN"],
        }
    )
    cfg = pl.PubChemCfg(delay=0)
    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] = {}
    parent_record_cache: dict[str, pd.Series | None] = {}

    (
        prepared_cache,
        prepared_resolution,
        prepared_parent_cache,
        pending_parent_ids,
        load_parent_record,
    ) = pipeline._prepare_pubchem_caches(
        frame,
        cfg,
        cache_path=None,
        cache_ttl_hours=None,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        parent_record_cache=parent_record_cache,
    )

    assert prepared_cache is cid_cache
    assert prepared_resolution is resolution_cache
    assert prepared_parent_cache is parent_record_cache
    assert pending_parent_ids == []
    assert "CHEMBL1" in prepared_parent_cache

    parent_series = prepared_parent_cache["CHEMBL1"]
    assert isinstance(parent_series, pd.Series)
    assert parent_series["canonical_smiles"] == "C"
    loaded_parent = load_parent_record("CHEMBL1")
    assert isinstance(loaded_parent, pd.Series)
    assert loaded_parent["canonical_smiles"] == "C"


def test_prefetch_parents_updates_cache(monkeypatch: pytest.MonkeyPatch) -> None:
    cfg = pl.PubChemCfg(delay=0, batch_size=2)
    parent_record_cache: dict[str, pd.Series | None] = {}
    fetched = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL2", "CHEMBL2"],
            "canonical_smiles": ["CC", "CC"],
        }
    )

    calls: list[Sequence[str]] = []

    def fake_get_testitem(*args: object, **kwargs: object) -> pd.DataFrame:
        chembl_ids = args[0] if args else kwargs.get("chembl_ids")
        calls.append(list(chembl_ids))
        return fetched

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)

    pipeline._prefetch_parents(
        ["CHEMBL2"],
        client=object(),
        api_cfg=ApiCfg(),
        cfg=cfg,
        timeout=5.0,
        testitem_fields=["canonical_smiles"],
        request_limit=3,
        parent_record_cache=parent_record_cache,
    )

    assert calls == [["CHEMBL2"]]
    assert "CHEMBL2" in parent_record_cache
    assert parent_record_cache["CHEMBL2"]["canonical_smiles"] == "CC"


def test_resolve_pubchem_cids_marks_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "canonical_smiles": ["C", "N", "O"],
        }
    )
    cfg = pl.PubChemCfg(delay=0, write_not_found_literal=True)
    cid_cache: dict[str, str | None] = {"CHEMBL2": "222"}
    resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] = {}
    calls: list[str] = []

    def fake_resolve(
        row: pd.Series,
        cache: MutableMapping[str, str | None],
        cfg: pl.PubChemCfg,
        **_: object,
    ) -> str | None:
        chembl_id = row["molecule_chembl_id"]
        calls.append(chembl_id)
        if chembl_id == "CHEMBL3":
            cache[chembl_id] = "333"
            return "333"
        cache[chembl_id] = None
        return None

    monkeypatch.setattr(pipeline, "resolve_pubchem_cid", fake_resolve)

    skip_mask = pd.Series(False, index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    cid_series, lookup_cids, cache_dirty = pipeline._resolve_pubchem_cids(
        frame,
        cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        load_parent_record=lambda _: None,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert list(cid_series) == ["Not Found", "222", "333"]
    assert lookup_cids == {"222", "333"}
    assert cache_dirty
    assert calls == ["CHEMBL1", "CHEMBL3"]


def test_resolve_pubchem_cids_skips_marked_rows(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "canonical_smiles": ["C", "N"],
        }
    )
    cfg = pl.PubChemCfg(delay=0)
    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] = {}
    calls: list[str] = []

    def fake_resolve(
        row: pd.Series,
        cache: MutableMapping[str, str | None],
        cfg: pl.PubChemCfg,
        **_: object,
    ) -> str:
        chembl_id = row["molecule_chembl_id"]
        calls.append(chembl_id)
        cache[chembl_id] = "111"
        return "111"

    monkeypatch.setattr(pipeline, "resolve_pubchem_cid", fake_resolve)

    skip_mask = pd.Series([True, False], index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    cid_series, lookup_cids, cache_dirty = pipeline._resolve_pubchem_cids(
        frame,
        cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        load_parent_record=lambda _: None,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert calls == ["CHEMBL2"]
    assert pd.isna(cid_series.iloc[0])
    assert cid_series.iloc[1] == "111"
    assert lookup_cids == {"111"}
    assert cache_dirty


def test_resolve_pubchem_cids_skips_cached_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "canonical_smiles": ["C", "N"],
        }
    )
    cfg = pl.PubChemCfg(delay=0)
    cid_cache: dict[str, str | None] = {"CHEMBL1": None}
    resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] = {}
    calls: list[str] = []

    def fake_resolve(
        row: pd.Series,
        cache: MutableMapping[str, str | None],
        cfg: pl.PubChemCfg,
        **_: object,
    ) -> str:
        chembl_id = row["molecule_chembl_id"]
        calls.append(chembl_id)
        cache[chembl_id] = "111"
        return "111"

    monkeypatch.setattr(pipeline, "resolve_pubchem_cid", fake_resolve)

    skip_mask = pd.Series(False, index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    cid_series, lookup_cids, cache_dirty = pipeline._resolve_pubchem_cids(
        frame,
        cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        load_parent_record=lambda _: None,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert calls == ["CHEMBL2"]
    assert pd.isna(cid_series.iloc[0])
    assert cid_series.iloc[1] == "111"
    assert lookup_cids == {"111"}
    assert cache_dirty


def test_merge_pubchem_properties_preserves_existing_values(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "pubchem_cid": ["111", "222", pd.NA],
            "pubchem_iupac_name": ["existing", "mixture", ""],
        }
    )
    cid_series = pd.Series(["111", "222", "333"], index=frame.index, dtype="string")
    lookup_cids = {"333"}
    cfg = pl.PubChemCfg(delay=0, prefer_local_values=True)

    props = {
        "333": pl.Properties("new", "formula", "iso", "can", "inchi", "inchikey"),
    }

    monkeypatch.setattr(pl, "get_properties", lambda cid, cfg: props[cid])

    skip_mask = pd.Series([True, True, False], index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    pubchem_df = pipeline._merge_pubchem_properties(
        frame,
        cid_series,
        lookup_cids,
        cfg=cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert pubchem_df.loc[0, "pubchem_iupac_name"] == "existing"
    assert pubchem_df.loc[1, "pubchem_iupac_name"] == "mixture"
    assert pubchem_df.loc[2, "pubchem_iupac_name"] == "new"
    assert pubchem_df.loc[2, "pubchem_molecular_formula"] == "formula"


def test_merge_pubchem_properties_throttles_by_rps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    count = 5
    frame = pd.DataFrame({"value": range(count)})
    cid_series = pd.Series(
        [f"CID-{idx}" for idx in range(count)], index=frame.index, dtype="string"
    )
    lookup_cids = set(cid_series.tolist())
    cfg = pl.PubChemCfg(delay=0, rps=2, batch_size=10)

    lock = threading.Lock()
    active = 0
    peak_active = 0

    def fake_get_properties(cid: str, pubchem_cfg: pl.PubChemCfg) -> pl.Properties:
        nonlocal active, peak_active
        assert pubchem_cfg is cfg
        with lock:
            active += 1
            peak_active = max(peak_active, active)
        time.sleep(0.05)
        with lock:
            active -= 1
        return pl.Properties(f"name-{cid}", None, None, None, None, None)

    monkeypatch.setattr(pl, "get_properties", fake_get_properties)

    skip_mask = pd.Series(False, index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    result = pipeline._merge_pubchem_properties(
        frame,
        cid_series,
        lookup_cids,
        cfg=cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert peak_active == cfg.rps
    assert result["pubchem_cid"].tolist() == cid_series.tolist()
    assert result["pubchem_iupac_name"].tolist() == [
        f"name-{cid}" for cid in cid_series.tolist()
    ]


def test_merge_pubchem_properties_limits_workers_to_rps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"value": [0, 1, 2]})
    cid_series = pd.Series(
        ["CID-0", "CID-1", "CID-2"], index=frame.index, dtype="string"
    )
    lookup_cids = set(cid_series.tolist())
    cfg = pl.PubChemCfg(delay=0, rps=2, batch_size=5)

    monkeypatch.setattr(
        pl,
        "get_properties",
        lambda cid, pubchem_cfg: pl.Properties(
            f"name-{cid}", None, None, None, None, None
        ),
    )

    class ImmediateFuture:
        def __init__(self, value: object) -> None:
            self._value = value

        def result(self) -> object:
            return self._value

    captured_workers: list[int] = []

    class DummyExecutor:
        def __init__(self, max_workers: int) -> None:
            captured_workers.append(max_workers)

        def __enter__(self) -> DummyExecutor:
            return self

        def __exit__(self, exc_type: object, exc: object, tb: object) -> bool:
            return False

        def submit(
            self, func: Callable[..., object], *args: object, **kwargs: object
        ) -> ImmediateFuture:
            return ImmediateFuture(func(*args, **kwargs))

    monkeypatch.setattr(pipeline, "ThreadPoolExecutor", DummyExecutor)

    skip_mask = pd.Series(False, index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    result = pipeline._merge_pubchem_properties(
        frame,
        cid_series,
        lookup_cids,
        cfg=cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert captured_workers == [cfg.rps]
    assert result["pubchem_iupac_name"].tolist() == [
        f"name-{cid}" for cid in cid_series.tolist()
    ]


def test_merge_pubchem_properties_uses_batch_size_when_below_rps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    count = 6
    frame = pd.DataFrame({"value": range(count)})
    cid_series = pd.Series(
        [f"CID-{idx}" for idx in range(count)], index=frame.index, dtype="string"
    )
    lookup_cids = set(cid_series.tolist())
    cfg = pl.PubChemCfg(delay=0, rps=10, batch_size=3)

    lock = threading.Lock()
    active = 0
    peak_active = 0

    def fake_get_properties(cid: str, pubchem_cfg: pl.PubChemCfg) -> pl.Properties:
        nonlocal active, peak_active
        assert pubchem_cfg is cfg
        with lock:
            active += 1
            peak_active = max(peak_active, active)
        time.sleep(0.05)
        with lock:
            active -= 1
        return pl.Properties(f"name-{cid}", None, None, None, None, None)

    monkeypatch.setattr(pl, "get_properties", fake_get_properties)

    skip_mask = pd.Series(False, index=frame.index)
    prefer_local_mask = pd.Series(False, index=frame.index)

    result = pipeline._merge_pubchem_properties(
        frame,
        cid_series,
        lookup_cids,
        cfg=cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert peak_active == cfg.batch_size
    assert result["pubchem_cid"].tolist() == cid_series.tolist()
    assert result["pubchem_iupac_name"].tolist() == [
        f"name-{cid}" for cid in cid_series.tolist()
    ]


def test_add_pubchem_data_fetches_properties_in_parallel(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": [
                "CHEMBL1",
                "CHEMBL2",
                "CHEMBL3",
                "CHEMBL4",
            ],
            "canonical_smiles": ["C", "CC", "CCC", "CCCC"],
        }
    )
    cfg = pl.PubChemCfg(delay=0, rps=2)

    monkeypatch.setattr(pipeline, "_load_pubchem_cid_cache", lambda *_, **__: {})
    monkeypatch.setattr(pipeline, "_write_pubchem_cid_cache", lambda *_, **__: None)

    cid_map = {
        "C": "1",
        "CC": "2",
        "CCC": "3",
        "CCCC": "4",
    }

    def fake_resolve(
        row: pd.Series,
        cache: MutableMapping[str, str | None],
        cfg: pl.PubChemCfg,
        **_: object,
    ) -> str:
        cid = cid_map[row["canonical_smiles"]]
        chembl_id = row["molecule_chembl_id"]
        if chembl_id:
            cache[chembl_id] = cid
        return cid

    monkeypatch.setattr(pipeline, "resolve_pubchem_cid", fake_resolve)

    lock = threading.Lock()
    active = 0
    peak_active = 0
    call_order: list[str] = []
    sleep_duration = 0.1

    def fake_get_properties(cid: str, pubchem_cfg: pl.PubChemCfg) -> pl.Properties:
        nonlocal active, peak_active
        with lock:
            active += 1
            peak_active = max(peak_active, active)
        call_order.append(cid)
        time.sleep(sleep_duration)
        with lock:
            active -= 1
        return pl.Properties(f"name-{cid}", None, None, None, None, None)

    monkeypatch.setattr(pl, "get_properties", fake_get_properties)

    start = time.perf_counter()
    result = pipeline.add_pubchem_data(df, cfg)
    elapsed = time.perf_counter() - start

    assert peak_active <= cfg.rps
    assert peak_active > 1
    assert sorted(call_order) == sorted(cid_map.values())
    assert elapsed < sleep_duration * len(cid_map) * 0.75

    expected_names = [f"name-{cid_map[value]}" for value in df["canonical_smiles"]]
    assert result["pubchem_iupac_name"].tolist() == expected_names


def test_resolve_pubchem_cid_prefers_inchikey(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    row = pd.Series(
        {
            "molecule_chembl_id": "chembl1",
            "standard_inchi_key": "abc-key",
            "standard_inchi": "ignored",
            "pref_name": "ignored",
            "canonical_smiles": "C",
        }
    )
    cfg = pl.PubChemCfg(delay=0)
    cache: dict[str, str | None] = {}
    captured: dict[str, str | None] = {}

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        cid_cache: Mapping[str, str | None] | None = None,
        cache_key: str | None,
        **kwargs: object,
    ) -> pl.PubChemResolution:
        captured.update(identifiers)
        assert cache_key == "CHEMBL1"
        if isinstance(cid_cache, dict) and cache_key:
            cid_cache[cache_key] = "10"
        return pl.PubChemResolution(cid="10", source="standard_inchi_key")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)

    cid = pipeline.resolve_pubchem_cid(row, cache, cfg, resolution_cache={})

    assert cid == "10"
    assert cache["CHEMBL1"] == "10"
    assert captured["standard_inchi_key"] == "ABC-KEY"


def test_resolve_pubchem_cid_uses_pubchem_inchikey(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    row = pd.Series(
        {
            "molecule_chembl_id": "chembl1",
            "standard_inchi_key": pd.NA,
            "pubchem_inchikey": "xyz-key",
            "pref_name": "ignored",
            "canonical_smiles": "C",
        }
    )
    cfg = pl.PubChemCfg(delay=0)
    cache: dict[str, str | None] = {}
    captured: dict[str, str | None] = {}

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        cid_cache: Mapping[str, str | None] | None = None,
        cache_key: str | None,
        **kwargs: object,
    ) -> pl.PubChemResolution:
        captured.update(identifiers)
        assert cache_key == "CHEMBL1"
        if isinstance(cid_cache, dict) and cache_key:
            cid_cache[cache_key] = "77"
        return pl.PubChemResolution(cid="77", source="pubchem_inchikey")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)

    cid = pipeline.resolve_pubchem_cid(row, cache, cfg, resolution_cache={})

    assert cid == "77"
    assert cache["CHEMBL1"] == "77"
    assert captured["pubchem_inchikey"] == "XYZ-KEY"


def test_resolve_pubchem_cid_uses_parent_when_enabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    child_row = pd.Series(
        {
            "molecule_chembl_id": "CHEMBL2",
            "parent_molecule_chembl_id": "CHEMBL1",
            "standard_inchi_key": pd.NA,
            "standard_inchi": pd.NA,
            "pref_name": pd.NA,
            "canonical_smiles": pd.NA,
        }
    )
    parent_row = pd.Series(
        {
            "molecule_chembl_id": "CHEMBL1",
            "standard_inchi_key": "parent-key",
        }
    )
    cache: dict[str, str | None] = {}
    cfg = pl.PubChemCfg(delay=0, use_parent_for_salts=True)

    calls: list[str] = []

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        cid_cache: Mapping[str, str | None] | None = None,
        cache_key: str | None,
        **kwargs: object,
    ) -> pl.PubChemResolution:
        calls.append(cache_key or "")
        if cache_key == "CHEMBL1":
            return pl.PubChemResolution(cid="42", source="standard_inchi_key")
        return pl.PubChemResolution(cid=None, source=None)

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)

    cid = pipeline.resolve_pubchem_cid(
        child_row,
        cache,
        cfg,
        parent_loader=lambda _: parent_row,
        resolution_cache={},
    )

    assert cid == "42"
    assert cache["CHEMBL1"] == "42"
    assert cache["CHEMBL2"] == "42"
    assert calls == ["CHEMBL2", "CHEMBL1"]


def test_resolve_pubchem_cid_handles_missing_parent_column(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    row = pd.Series({"molecule_chembl_id": "CHEMBL2"})
    cache: dict[str, str | None] = {}
    cfg = pl.PubChemCfg(delay=0, use_parent_for_salts=True)

    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: pl.PubChemResolution(cid=None, source=None),
    )

    cid = pipeline.resolve_pubchem_cid(row, cache, cfg, resolution_cache={})

    assert cid is None
    assert cache["CHEMBL2"] is None


def test_resolve_pubchem_cid_logs_when_parent_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    row = pd.Series(
        {
            "molecule_chembl_id": "CHEMBL2",
            "parent_molecule_chembl_id": "CHEMBL1",
        }
    )
    cache: dict[str, str | None] = {}
    cfg = pl.PubChemCfg(delay=0, use_parent_for_salts=True)

    events: list[tuple[str, dict[str, object]]] = []

    def capture(event: str, **kwargs: object) -> None:
        events.append((event, kwargs))

    monkeypatch.setattr(pipeline.logger, "info", capture)
    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: pl.PubChemResolution(cid=None, source=None),
    )

    cid = pipeline.resolve_pubchem_cid(
        row,
        cache,
        cfg,
        parent_loader=lambda _: None,
        resolution_cache={},
    )

    assert cid is None
    assert cache["CHEMBL2"] is None
    assert any(event == "pubchem_parent_structure_missing" for event, _ in events)


def test_prepare_pubchem_caches_uses_disk_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_path = tmp_path / "cid_cache.json"
    cache_path.write_text(json.dumps({"CHEMBL1": "321"}))

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    cfg = pl.PubChemCfg(delay=0, cid_cache_path=cache_path)

    (
        cid_cache,
        resolution_cache,
        parent_record_cache,
        pending_parent_ids,
        load_parent_record,
    ) = pipeline._prepare_pubchem_caches(
        frame,
        cfg,
        cache_path=cache_path,
        cache_ttl_hours=None,
        cid_cache=None,
        resolution_cache=None,
        parent_record_cache=None,
    )

    assert cid_cache.get("CHEMBL1") == "321"
    assert isinstance(resolution_cache, dict)
    assert isinstance(parent_record_cache, dict)
    assert pending_parent_ids == []
    assert load_parent_record("CHEMBL1") is not None


def test_write_pubchem_cid_cache_creates_parent_dir(tmp_path: Path) -> None:
    cache_path = tmp_path / "nested" / "cid_cache.json"

    assert not cache_path.parent.exists()

    pipeline._write_pubchem_cid_cache(cache_path, {"CHEMBL1": "321", "CHEMBL2": None})

    assert cache_path.exists()
    payload = json.loads(cache_path.read_text())
    assert payload["values"] == {"CHEMBL1": "321", "CHEMBL2": None}
    assert payload["metadata"]["version"] == pipeline._PUBCHEM_CACHE_SCHEMA_VERSION


def test_write_pubchem_cid_cache_partial_write_keeps_original(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_path = tmp_path / "cid_cache.json"
    original_payload = {
        "metadata": {
            "version": pipeline._PUBCHEM_CACHE_SCHEMA_VERSION,
            "updated_at": datetime.now(UTC).isoformat(),
        },
        "values": {"CHEMBL1": "321"},
    }
    cache_path.write_text(json.dumps(original_payload))

    real_dump = pipeline.json.dump

    def fake_dump(
        payload: object, handle: Any, *args: object, **kwargs: object
    ) -> None:
        if isinstance(payload, dict) and "values" in payload:
            handle.write('{"values": {')
            handle.flush()
            raise OSError("disk full")
        real_dump(payload, handle, *args, **kwargs)

    monkeypatch.setattr(pipeline.json, "dump", fake_dump)

    pipeline._write_pubchem_cid_cache(cache_path, {"CHEMBL2": "654"})

    assert json.loads(cache_path.read_text()) == original_payload


def test_load_pubchem_cid_cache_uses_shared_selector(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_path = tmp_path / "cid_cache.json"
    cache_payload = {
        "values": {
            "chembl1": "123 | 456",
            "chembl2": None,
        }
    }
    cache_path.write_text(json.dumps(cache_payload))

    calls: list[tuple[str | None, str, str | None, str | None]] = []

    def fake_select(
        candidates: str | None,
        *,
        identifier: str,
        value: str | None,
        context_id: str | None = None,
    ) -> str | None:
        calls.append((candidates, identifier, value, context_id))
        if candidates:
            return candidates.split("|")[0].strip()
        return None

    monkeypatch.setattr(pl, "select_primary_cid", fake_select)

    cache = gtd._load_pubchem_cid_cache(cache_path)

    assert cache == {"CHEMBL1": "123", "CHEMBL2": None}
    assert calls == [("123 | 456", "cache_file", "123 | 456", "CHEMBL1")]


def test_load_pubchem_cid_cache_expires(tmp_path: Path) -> None:
    cache_path = tmp_path / "cid_cache.json"
    expired_at = datetime.now(UTC) - timedelta(hours=2)
    payload = {
        "metadata": {
            "version": pipeline._PUBCHEM_CACHE_SCHEMA_VERSION,
            "updated_at": expired_at.isoformat(),
        },
        "values": {"CHEMBL1": "321"},
    }
    cache_path.write_text(json.dumps(payload))

    cache = pipeline._load_pubchem_cid_cache(cache_path, ttl_hours=1)

    assert cache == {}


def test_run_chembl_column_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    """Ensure schema columns precede alphabetically sorted extras."""
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL1",
                "parent_molecule_id": "CHEMBL0",
                "molecule_type": "Small molecule",
                "salt_chembl_id": "CHEMBL1-SALT",
                "chirality": 1,
                "extra_b": 2,
                "extra_a": 1,
            }
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df)
    precomputed_catalog = {"CHEMBL1": "CHEMBL1_PARENT"}
    monkeypatch.setattr(
        pipeline, "load_parent_catalog", lambda **__: precomputed_catalog
    )
    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})

    monkeypatch.setattr(
        pipeline,
        "add_pubchem_data",
        lambda df, cfg, **__: df,
    )
    captured_precomputed: dict[str, object] = {"data": None, "frame": None}

    def fake_attach_parent_molecule_ids(frame: pd.DataFrame, **kwargs: object):
        captured_precomputed["data"] = kwargs.get("precomputed")
        captured_precomputed["frame"] = frame.copy()
        return (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        )

    monkeypatch.setattr(
        pipeline, "attach_parent_molecule_ids", fake_attach_parent_molecule_ids
    )
    monkeypatch.setattr(
        pipeline, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda p: "deadbeef")

    captured: dict[str, list[str]] = {}

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        captured["col_order"] = list(col_order or [])
        captured["columns"] = list(df.columns)
        return output

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0

    prepared = captured_precomputed["data"]
    assert isinstance(prepared, pipeline.ParentLookupPreparedData)
    captured_frame = captured_precomputed["frame"]
    assert isinstance(captured_frame, pd.DataFrame)
    expected_prepared = prepare_parent_lookup_data(captured_frame, cfg.molecule_catalog)
    pd.testing.assert_series_equal(
        prepared.child_ids, expected_prepared.child_ids, check_names=False
    )
    pd.testing.assert_series_equal(
        prepared.existing_parent_ids,
        expected_prepared.existing_parent_ids,
        check_names=False,
    )
    assert (
        prepared.existing_parent_ids.tolist()
        == expected_prepared.existing_parent_ids.tolist()
    )
    parent_field = cfg.molecule_catalog.parent_field
    child_field = cfg.molecule_catalog.child_field
    existing_parent = (
        captured_frame.get(parent_field)
        if parent_field in captured_frame
        else pd.Series([""], index=captured_frame.index, dtype="string")
    )
    normalised_parent = existing_parent.astype("string").fillna("").str.strip()
    child_ids = captured_frame.get(child_field, pd.Series([], dtype="string"))
    if child_ids.empty:
        child_ids = prepared.child_ids
    normalised_child = child_ids.astype("string").fillna("").str.upper()
    expected_need_lookup = set(normalised_child[normalised_parent == ""])
    expected_need_lookup -= set(precomputed_catalog)
    assert prepared.need_lookup == expected_need_lookup

    available = set(captured.get("columns", []))
    expected_head = [c for c in TestitemsSchema.columns if c in available]
    expected_tail = sorted(available - set(TestitemsSchema.columns))
    assert captured["col_order"] == expected_head + expected_tail


def test_run_chembl_parent_lookup_precomputed_excludes_resolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\nCHEMBL2\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1", "CHEMBL2"]))

    df = pd.DataFrame(
        [
            {"molecule_chembl_id": "CHEMBL1", "parent_molecule_id": pd.NA},
            {"molecule_chembl_id": "CHEMBL2", "parent_molecule_id": pd.NA},
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df.copy())
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})

    cache_path = tmp_path / "parent_catalog.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"

    parent_catalog = {"CHEMBL1": "CHEMBL1_PARENT"}

    monkeypatch.setattr(
        pipeline, "query_parent_catalog", lambda *_, **__: parent_catalog
    )
    monkeypatch.setattr(
        pipeline.molecule_catalog, "fetch_parent_catalog_for", lambda *_, **__: {}
    )

    captured_precomputed: dict[str, object] = {}

    def fake_attach(frame: pd.DataFrame, **kwargs: object):
        captured_precomputed["data"] = kwargs.get("precomputed")
        return (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
                missing=0,
                unique=0,
                attached=len(frame),
                uncovered=0,
            ),
        )

    monkeypatch.setattr(pipeline, "attach_parent_molecule_ids", fake_attach)
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0

    prepared = captured_precomputed.get("data")
    assert isinstance(prepared, pipeline.ParentLookupPreparedData)
    assert prepared.need_lookup == {"CHEMBL2"}


def test_run_chembl_initialises_pubchem_session(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        {"molecule_chembl_id": ["CHEMBL1"], "molecule_type": ["Small molecule"]}
    )
    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df)
    monkeypatch.setattr(
        pipeline, "add_pubchem_data", lambda frame, pubchem_cfg, **__: frame
    )
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})

    monkeypatch.setattr(
        pipeline,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )

    captured: dict[str, object] = {}

    def fake_init_session(api: object, retry: object) -> None:
        captured["init"] = (api, retry)

    monkeypatch.setattr(gtd.pc, "init_session", fake_init_session)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(
        pipeline, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0

    init_api, init_retry = captured["init"]
    expected_api = gtd._prepare_pubchem_api_cfg(cfg, cfg.api)

    assert init_api == expected_api
    assert init_retry == cfg.retry


def test_run_chembl_uses_lazy_identifier_stream(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text(
        "molecule_chembl_id\nCHEMBL1\nCHEMBL2\n",
        encoding=cfg.io.csv_encoding,
    )

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    class SingleUseIterable:
        def __init__(self, values: list[str]) -> None:
            self._values = values
            self.iterations = 0

        def __iter__(self):
            if self.iterations:
                raise AssertionError("identifiers iterator consumed multiple times")
            self.iterations += 1
            yield from self._values

    ids_source = SingleUseIterable(["CHEMBL1", "CHEMBL2"])

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: ids_source)

    received_ids: object | None = None

    def fake_get_testitem(ids, **kwargs):  # type: ignore[no-untyped-def]
        nonlocal received_ids
        received_ids = ids
        values = list(ids)
        return pd.DataFrame(
            [
                {
                    cfg.molecule_catalog.child_field: value,
                    cfg.molecule_catalog.parent_field: pd.NA,
                }
                for value in values
            ]
        )

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(gtd.pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, *args, **__: frame)
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(pipeline, "query_parent_catalog", lambda *_, **__: {})
    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})
    monkeypatch.setattr(
        pipeline,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )
    monkeypatch.setattr(pipeline, "normalize_testitems", lambda frame: frame)
    monkeypatch.setattr(pipeline, "add_pipeline_metadata", lambda frame: frame)
    monkeypatch.setattr(
        pipeline,
        "validate_testitems",
        lambda frame, return_result=True: SimpleNamespace(
            data=frame,
            failure_cases=pd.DataFrame(),
        ),
    )
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "update_parent_catalog_cache", lambda *_, **__: None)

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert ids_source.iterations == 1
    assert received_ids is not None
    assert not isinstance(received_ids, list)
    assert isinstance(received_ids, Iterator)
    assert received_ids.__class__.__name__ == "generator"


def test_run_chembl_calls_pubchem_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    df = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL1",
                "canonical_smiles": "C",
                "molecule_type": "Small molecule",
                "parent_molecule_chembl_id": pd.NA,
            }
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: df.copy())
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})

    monkeypatch.setattr(
        pipeline,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )

    cfg.pubchem.cid_cache_path = tmp_path / "pubchem_cache.json"
    cfg.pubchem.delay = 0

    resolve_calls: list[Mapping[str, str | None]] = []

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        pubchem_cfg: pl.PubChemCfg,
        **kwargs: object,
    ) -> pl.PubChemResolution:
        resolve_calls.append(identifiers)
        return pl.PubChemResolution(cid="321", source="resolver")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)

    props = pl.Properties("name", "formula", "i", "c", "inchi", "inchikey")
    monkeypatch.setattr(pl, "get_properties", lambda cid, cfg: props)
    monkeypatch.setattr(gtd.pc, "init_session", lambda api, retry: None)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda frame, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(
        pipeline, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert len(resolve_calls) == 1


def test_run_chembl_prefills_parent_from_hierarchy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text(
        "molecule_chembl_id\nCHEMBL1\nCHEMBL2\n", encoding=cfg.io.csv_encoding
    )

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1", "CHEMBL2"]))

    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    source = pd.DataFrame(
        [
            {child_field: "CHEMBL1", parent_field: pd.NA},
            {child_field: "CHEMBL2", parent_field: "CHEMBL2_EXISTING"},
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: source.copy())
    monkeypatch.setattr(gtd.pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **_: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(pipeline, "update_parent_catalog_cache", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_parent_catalog_cache", lambda *_, **__: None)

    cfg.molecule_catalog.cache_path = tmp_path / "parent_catalog.json"
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"
    hierarchy_path = tmp_path / "hierarchy.csv"
    cfg.testitem_molecule_enrichment.sources.molecule_hierarchy_path = hierarchy_path

    captured_path: Path | None = None

    def fake_lookup(
        path: Path | None, *, io_cfg: object, **_: object
    ) -> dict[str, str]:
        nonlocal captured_path
        captured_path = path
        return {"CHEMBL1": "CHEMBL1_PARENT"}

    def fail_query(*_: object, **__: object) -> dict[str, str]:
        raise AssertionError("query_parent_catalog should not be called")

    def fail_fetch(*_: object, **__: object) -> dict[str, str]:
        raise AssertionError("fetch_parent_catalog_for should not be called")

    def fail_load(**_: object) -> dict[str, str]:
        raise AssertionError("load_parent_catalog should not be called")

    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", fake_lookup)
    monkeypatch.setattr(pipeline, "query_parent_catalog", fail_query)
    monkeypatch.setattr(
        pipeline.molecule_catalog, "fetch_parent_catalog_for", fail_fetch
    )
    monkeypatch.setattr(pipeline, "load_parent_catalog", fail_load)

    captured_df: pd.DataFrame | None = None

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        nonlocal captured_df
        captured_df = df.copy()
        return output

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert captured_path == hierarchy_path
    assert captured_df is not None
    assert captured_df[parent_field].tolist() == ["CHEMBL1_PARENT", "CHEMBL2_EXISTING"]


def test_run_chembl_merges_parent_catalog(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\nCHEMBL2\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1", "CHEMBL2"]))

    source = pd.DataFrame(
        [
            {"molecule_chembl_id": "CHEMBL1", "parent_molecule_chembl_id": None},
            {
                "molecule_chembl_id": "CHEMBL2",
                "parent_molecule_chembl_id": "CHEMBL2_EXISTING",
            },
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: source.copy())
    cache_path = tmp_path / "parent_catalog.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"

    parent_catalog = {"CHEMBL1": "CHEMBL1_PARENT", "CHEMBL2": "CHEMBL2_PARENT"}
    query_calls = 0

    def fake_query_parent_catalog(
        ids: Iterable[str], *, catalog_cfg: object
    ) -> dict[str, str]:
        nonlocal query_calls
        query_calls += 1
        return parent_catalog

    monkeypatch.setattr(pipeline, "query_parent_catalog", fake_query_parent_catalog)
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})
    captured_catalog: dict[str, str] | None = None
    captured_source: str | None = None

    def fake_attach_parent_molecule_ids(
        frame: pd.DataFrame, **kwargs: object
    ) -> tuple[pd.DataFrame, pipeline.ParentLookupStats]:
        nonlocal captured_catalog, captured_source
        captured_catalog = kwargs.get("catalog")
        captured_source = kwargs.get("source")

        catalog_cfg = kwargs.get("catalog_cfg")
        mapping = captured_catalog or {}
        parent_field = (
            catalog_cfg.parent_field
            if catalog_cfg is not None
            else "parent_molecule_chembl_id"
        )
        child_field = (
            catalog_cfg.child_field if catalog_cfg is not None else "molecule_chembl_id"
        )
        updated = frame.copy()
        parent_series = updated[parent_field].astype("string")
        mask = parent_series.isna() | parent_series.eq("")
        updated.loc[mask, parent_field] = (
            updated.loc[mask, child_field].map(mapping).astype("string")
        )
        return (
            updated,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        )

    monkeypatch.setattr(
        pipeline, "attach_parent_molecule_ids", fake_attach_parent_molecule_ids
    )
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    captured_df: pd.DataFrame | None = None

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        nonlocal captured_df
        captured_df = df.copy()
        return output

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0
    assert captured_df is not None
    assert captured_df["parent_molecule_chembl_id"].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_EXISTING",
    ]
    assert query_calls == 1
    assert captured_catalog is parent_catalog
    assert captured_source == pipeline.PARENT_LOOKUP_SOURCE_LOOKUP


def test_run_chembl_updates_parent_cache_and_reuses_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding)

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    child_field = cfg.molecule_catalog.child_field
    parent_field = cfg.molecule_catalog.parent_field
    source = pd.DataFrame([{child_field: "CHEMBL1", parent_field: pd.NA}])

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: source.copy())
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})

    cache_path = tmp_path / "parent_catalog.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"

    fetch_calls: list[list[str]] = []

    def fake_fetch(
        ids: Iterable[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        normalised = sorted(ids)
        fetch_calls.append(normalised)
        if len(fetch_calls) > 1:
            raise AssertionError("fetch_parent_catalog_for should not be called twice")
        return {item: f"{item}_PARENT" for item in normalised}

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    update_calls: list[dict[str, str]] = []
    original_update = pipeline.update_parent_catalog_cache

    def tracking_update(catalog: Mapping[str, str], catalog_cfg: object) -> None:
        update_calls.append(dict(catalog))
        original_update(catalog, catalog_cfg)

    monkeypatch.setattr(pipeline, "update_parent_catalog_cache", tracking_update)

    attach_calls: list[tuple[dict[str, str], str | None]] = []

    def fake_attach(
        frame: pd.DataFrame, **kwargs: object
    ) -> tuple[pd.DataFrame, pipeline.ParentLookupStats]:
        catalog = dict(kwargs.get("catalog", {}))
        source_name = kwargs.get("source")
        attach_calls.append((catalog, source_name))
        stats_source = source_name or pipeline.PARENT_LOOKUP_SOURCE_SKIPPED
        stats = pipeline.ParentLookupStats(
            source=stats_source,
            missing=0,
            unique=len(catalog),
            attached=len(frame),
            uncovered=0,
        )
        return frame, stats

    monkeypatch.setattr(pipeline, "attach_parent_molecule_ids", fake_attach)
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )

    rc_first = gtd.run_chembl(cfg, args)
    rc_second = gtd.run_chembl(cfg, args)

    assert rc_first == 0
    assert rc_second == 0
    assert fetch_calls == [["CHEMBL1"]]
    assert update_calls == [{"CHEMBL1": "CHEMBL1_PARENT"}]
    assert [source for _, source in attach_calls] == [
        pipeline.PARENT_LOOKUP_SOURCE_PARTIAL,
        pipeline.PARENT_LOOKUP_SOURCE_LOOKUP,
    ]


def test_attach_parent_ids_preserves_existing_values_when_cache_has_no_matches(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
) -> None:
    parent_field = cfg.molecule_catalog.parent_field
    cache_path = tmp_path / "parent_cache.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path

    source = pd.DataFrame(
        {
            cfg.molecule_catalog.child_field: ["CHEMBL1", "CHEMBL2"],
            parent_field: ["CHEMBL1_EXISTING", pd.NA],
        }
    )

    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    result, stats = pipeline.attach_parent_molecule_ids(
        source,
        client=object(),
        api_cfg=cfg.api,
        catalog_cfg=cfg.molecule_catalog,
        timeout=None,
    )

    expected = pd.Series(
        ["CHEMBL1_EXISTING", pd.NA], index=result.index, dtype="string"
    )
    pd.testing.assert_series_equal(result[parent_field], expected, check_names=False)
    assert stats.missing == 1
    assert stats.attached == 1


def test_run_chembl_preserves_existing_parent_value_when_catalog_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    source = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL1",
                cfg.molecule_catalog.parent_field: "CHEMBL1_EXISTING",
            }
        ]
    )

    monkeypatch.setattr(cl, "get_testitem", lambda *_, **__: source.copy())
    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})
    monkeypatch.setattr(
        pipeline,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    captured_df: pd.DataFrame | None = None

    def fake_write_csv(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config,
        key_cols: list[str] | None = None,
        col_order: list[str] | None = None,
        **__: object,
    ) -> Path:
        nonlocal captured_df
        captured_df = df.copy()
        return output

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0
    assert captured_df is not None
    assert captured_df[cfg.molecule_catalog.parent_field].tolist() == [
        "CHEMBL1_EXISTING",
    ]


def test_run_chembl_parent_catalog_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    cache_path = tmp_path / "parent_catalog.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"

    monkeypatch.setattr(
        cl,
        "get_testitem",
        lambda *_, **__: pd.DataFrame(
            [
                {
                    cfg.molecule_catalog.child_field: "CHEMBL1",
                    cfg.molecule_catalog.parent_field: pd.NA,
                }
            ]
        ),
    )
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        pipeline,
        "query_parent_catalog",
        lambda *_, **__: (_ for _ in ()).throw(
            ValueError("missing columns: parant_molecule_id")
        ),
    )

    called = False

    def fake_write_csv(*args: object, **kwargs: object) -> Path:  # pragma: no cover
        nonlocal called
        called = True
        return Path("unused.csv")

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 1
    assert not called


def test_run_chembl_parent_catalog_request_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text("molecule_chembl_id\nCHEMBL1\n")

    args = argparse.Namespace(input_csv=input_csv, output_csv=tmp_path / "out.csv")

    monkeypatch.setattr(io, "read_ids", lambda *_, **__: iter(["CHEMBL1"]))

    cache_path = tmp_path / "parent_catalog.json"
    cache_path.write_text("{}", encoding="utf-8")
    cfg.molecule_catalog.cache_path = cache_path
    cfg.molecule_catalog.sqlite_path = tmp_path / "parent_catalog.sqlite"

    monkeypatch.setattr(
        cl,
        "get_testitem",
        lambda *_, **__: pd.DataFrame(
            [
                {
                    cfg.molecule_catalog.child_field: "CHEMBL1",
                    cfg.molecule_catalog.parent_field: pd.NA,
                }
            ]
        ),
    )
    monkeypatch.setattr(pipeline, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(pipeline, "load_molecule_hierarchy_lookup", lambda *_, **__: {})
    monkeypatch.setattr(pipeline, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(pipeline, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        pipeline,
        "query_parent_catalog",
        lambda *_, **__: (_ for _ in ()).throw(requests.RequestException("boom")),
    )

    monkeypatch.setattr(
        pipeline,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )

    called = False

    def fake_write_csv(*args: object, **kwargs: object) -> Path:  # pragma: no cover
        nonlocal called
        called = True
        return Path("unused.csv")

    monkeypatch.setattr(pipeline, "write_csv_deterministic", fake_write_csv)

    errors: list[tuple[str, dict[str, object]]] = []

    def fake_logger_error(event: str, *args: object, **kwargs: object) -> None:
        errors.append((event, kwargs))

    monkeypatch.setattr(pipeline.logger, "error", fake_logger_error)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 1
    assert not called
    assert any(
        event == "parent_catalog_invalid"
        and details.get("error") == "boom"
        and details.get("path") == str(cfg.molecule_catalog.cache_path)
        for event, details in errors
    )


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_fetches_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    monkeypatch.setattr(
        pipeline,
        "load_parent_catalog",
        lambda **__: {"CHEMBL1": "CHEMBL1_PARENT"},
    )

    fetched: dict[str, str] = {"CHEMBL2": "CHEMBL2_PARENT"}
    captured_ids: dict[str, list[str]] = {}

    def fake_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        captured_ids["ids"] = list(ids)
        return fetched

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    assert captured_ids["ids"] == ["CHEMBL1", "CHEMBL2"]
    assert result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert stats.unique == 2
    assert stats.attached == 2
    assert stats.missing == 0
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_SYNC


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_skips_filled_children(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    child_field = cfg.sources.chembl.molecule_catalog.child_field
    parent_field = cfg.sources.chembl.molecule_catalog.parent_field
    df = pd.DataFrame(
        {
            child_field: ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            parent_field: ["CHEMBL1_PARENT", pd.NA, "CHEMBL3_PARENT"],
        }
    )

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    query_calls: list[tuple[str, ...]] = []
    fetch_calls: list[list[str]] = []

    def fake_query(children: Iterable[str], catalog_cfg: object) -> dict[str, str]:
        captured = tuple(str(child) for child in children)
        query_calls.append(captured)
        return {}

    def fake_fetch(
        ids: Iterable[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        ordered = [str(child) for child in ids]
        fetch_calls.append(ordered)
        return {value: f"{value}_PARENT" for value in ordered}

    monkeypatch.setattr(pipeline, "query_parent_catalog", fake_query)
    monkeypatch.setattr(
        pipeline.molecule_catalog, "fetch_parent_catalog_for", fake_fetch
    )
    monkeypatch.setattr(pipeline, "write_parent_catalog_cache", lambda *_, **__: None)
    monkeypatch.setattr(pipeline, "update_parent_catalog_cache", lambda *_, **__: None)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    expected_lookup = ("CHEMBL2",)
    assert query_calls == [expected_lookup]
    assert fetch_calls == [["CHEMBL2"]]
    assert result[parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
        "CHEMBL3_PARENT",
    ]
    assert stats.unique == 1
    assert stats.attached == 3
    assert stats.missing == 0
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_PARTIAL


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_prefers_partial_fetch_when_complete(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    child_field = cfg.sources.chembl.molecule_catalog.child_field
    parent_field = cfg.sources.chembl.molecule_catalog.parent_field
    df = pd.DataFrame({child_field: ["CHEMBL10", "CHEMBL11"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    fetch_calls: list[list[str]] = []

    def fake_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        ordered = sorted(ids)
        fetch_calls.append(ordered)
        return {value: f"{value}_PARENT" for value in ordered}

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    monkeypatch.setattr(
        pipeline,
        "load_parent_catalog",
        lambda **_: (_ for _ in ()).throw(
            AssertionError("load_parent_catalog should not be called")
        ),
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    assert fetch_calls == [["CHEMBL10", "CHEMBL11"]]
    assert result[parent_field].tolist() == ["CHEMBL10_PARENT", "CHEMBL11_PARENT"]
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_PARTIAL


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_handles_large_catalog(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        }
    )

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    base_catalog = {f"CHEMBL{i}": f"CHEMBL{i}_PARENT" for i in range(1, 5000)}

    class TrackingMapping(Mapping[str, str]):
        def __init__(self, data: Mapping[str, str]) -> None:
            self._data = dict(data)
            self.lookups: list[str] = []
            self.iterations = 0

        def __getitem__(self, key: str) -> str:
            self.lookups.append(key)
            return self._data[key]

        def __iter__(self):  # type: ignore[override]
            self.iterations += 1
            return iter(self._data)

        def __len__(self) -> int:
            return len(self._data)

        def __contains__(self, key: object) -> bool:
            return key in self._data

    tracking_catalog = TrackingMapping(base_catalog)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        catalog=tracking_catalog,
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        precomputed=precomputed,
    )

    expected_values = [f"CHEMBL{i}_PARENT" for i in range(1, 4)]
    assert result[catalog_cfg.parent_field].tolist() == expected_values
    assert stats.unique == 3
    assert stats.attached == 3
    assert stats.missing == 0
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_LOOKUP
    assert tracking_catalog.iterations == 0
    assert set(tracking_catalog.lookups) == {"CHEMBL1", "CHEMBL2", "CHEMBL3"}


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_handles_partial_remote_success(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    monkeypatch.setattr(pipeline, "load_parent_catalog", lambda **__: {})

    remote_result = {"CHEMBL1": "CHEMBL1_PARENT"}

    def fake_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        assert ids == ["CHEMBL1", "CHEMBL2"]
        return remote_result

    monkeypatch.setattr(
        pipeline.molecule_catalog, "fetch_parent_catalog_for", fake_fetch
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    parent_values = result[catalog_cfg.parent_field].tolist()
    assert parent_values == ["CHEMBL1_PARENT", pd.NA]
    assert stats.attached == 1
    assert stats.missing == 1
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_SYNC


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_updates_cache_for_reuse(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.cache_path.parent.mkdir(parents=True, exist_ok=True)
    catalog_cfg.cache_path.write_text(
        json.dumps({"CHEMBL1": "CHEMBL1_PARENT"}, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    fetch_calls: list[list[str]] = []

    def fake_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        fetch_calls.append(list(ids))
        return {"CHEMBL2": "CHEMBL2_PARENT"}

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    first_result, first_stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    assert fetch_calls == [["CHEMBL2"]]
    assert first_result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert first_stats.source == pipeline.PARENT_LOOKUP_SOURCE_PARTIAL

    stored_catalog = json.loads(catalog_cfg.cache_path.read_text(encoding="utf-8"))
    assert stored_catalog == {
        "CHEMBL1": "CHEMBL1_PARENT",
        "CHEMBL2": "CHEMBL2_PARENT",
    }

    precomputed_second = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    second_result, second_stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed_second,
    )

    assert fetch_calls == [["CHEMBL2"]]
    assert second_result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert second_stats.source == pipeline.PARENT_LOOKUP_SOURCE_CACHE


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_uses_sqlite_after_migration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"
    catalog_cfg.cache_path.write_text(
        json.dumps({"CHEMBL1": "CHEMBL1_PARENT"}, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    original_loads = pipeline.molecule_catalog.json.loads
    calls = {"loads": 0}

    def counting_loads(data: str, *args: object, **kwargs: object) -> object:
        calls["loads"] += 1
        return original_loads(data, *args, **kwargs)

    monkeypatch.setattr(pipeline.molecule_catalog.json, "loads", counting_loads)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    first_result, _ = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    assert calls["loads"] == 1
    assert first_result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]

    calls["loads"] = 0

    precomputed_second = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    second_result, _ = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed_second,
    )

    assert calls["loads"] == 0
    assert second_result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_fetch_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    monkeypatch.setattr(
        pipeline,
        "load_parent_catalog",
        lambda **__: {},
    )

    def failing_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
        catalog_cfg: object | None = None,
    ) -> dict[str, str]:
        raise requests.RequestException("boom")

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        failing_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    parent_values = result[catalog_cfg.parent_field].tolist()
    assert parent_values == [pd.NA]
    assert stats.unique == 1
    assert stats.attached == 0
    assert stats.missing == 1
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_SYNC


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_full_sync_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"

    def failing_load_parent_catalog(**_: object) -> dict[str, str]:
        raise requests.RequestException("boom")

    monkeypatch.setattr(pipeline, "load_parent_catalog", failing_load_parent_catalog)

    target_catalog = getattr(pipeline, "molecule_catalog", molecule_catalog)
    monkeypatch.setattr(
        target_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )

    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    parent_values = result[catalog_cfg.parent_field].tolist()
    assert parent_values == [pd.NA]
    assert stats.unique == 1
    assert stats.attached == 0
    assert stats.missing == 1
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_SYNC


@pytest.mark.parametrize("use_precomputed", [False, True])
def test_attach_parent_molecule_ids_uses_cache_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cfg: Config,
    use_precomputed: bool,
) -> None:
    df = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    catalog_cfg = cfg.sources.chembl.molecule_catalog.model_copy(deep=True)
    catalog_cfg.cache_path = tmp_path / "catalog.json"
    catalog_cfg.sqlite_path = tmp_path / "catalog.sqlite"
    catalog_cfg.cache_path.write_text("{}", encoding="utf-8")

    def unexpected_load_parent_catalog(**_: object) -> dict[str, str]:
        raise AssertionError("load_parent_catalog should not be called")

    monkeypatch.setattr(pipeline, "load_parent_catalog", unexpected_load_parent_catalog)

    def unexpected_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
    ) -> dict[str, str]:
        raise AssertionError("fetch_parent_catalog_for should not be called")

    monkeypatch.setattr(
        pipeline.molecule_catalog,
        "fetch_parent_catalog_for",
        unexpected_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = pipeline.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        catalog={"CHEMBL1": "CHEMBL1_PARENT"},
        source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
        precomputed=precomputed,
    )

    assert result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]
    assert stats.missing == 0
    assert stats.attached == 1
    assert stats.source == pipeline.PARENT_LOOKUP_SOURCE_LOOKUP
