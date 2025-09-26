from __future__ import annotations

import argparse
import json
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Callable, Iterable, Mapping
from types import SimpleNamespace

import pandas as pd
import pytest
import requests

from library import chembl_library as cl
from library import io
from library import pubchem_library as pl
from library.config import Config
from schemas import TestitemsSchema
from scripts import get_testitem_data as gtd


def prepare_parent_lookup_data(
    df: pd.DataFrame, catalog_cfg
) -> gtd.ParentLookupPreparedData:
    child_column = catalog_cfg.child_field
    parent_column = catalog_cfg.parent_field

    if child_column in df.columns:
        normalised_child = gtd._normalise_chembl_ids(df[child_column])
    else:
        normalised_child = pd.Series("", index=df.index, dtype="string")

    if parent_column in df.columns:
        existing_parent = gtd._normalise_chembl_ids(df[parent_column])
    else:
        existing_parent = pd.Series("", index=df.index, dtype="string")

    need_lookup_mask = (normalised_child != "") & (existing_parent == "")
    need_lookup = set(normalised_child[need_lookup_mask])

    return gtd.ParentLookupPreparedData(
        child_ids=normalised_child,
        existing_parent_ids=existing_parent,
        need_lookup=need_lookup,
    )


def test_add_pubchem_data_missing_uses_na(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"canonical_smiles": ["C"]})
    cfg = pl.PubChemCfg(delay=0)

    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: pl.PubChemResolution(cid=None, source=None),
    )

    result = gtd.add_pubchem_data(df, cfg)
    pubchem_cols = [col for col in result.columns if col.startswith("pubchem_")]
    assert pubchem_cols
    assert result[pubchem_cols].isna().all().all()
    assert not (result[pubchem_cols] == "Not Found").any().any()


def test_add_pubchem_data_disabled(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"canonical_smiles": ["C"]})
    cfg = pl.PubChemCfg(delay=0, enable=False)

    calls: list[tuple[tuple, dict]] = []
    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: calls.append((args, kwargs)),
    )

    result = gtd.add_pubchem_data(df, cfg)

    assert calls == []
    assert list(result.columns) == list(df.columns)


def test_add_pubchem_data_not_found_literal(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"canonical_smiles": ["C"]})
    cfg = pl.PubChemCfg(delay=0, write_not_found_literal=True)

    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: pl.PubChemResolution(cid=None, source=None),
    )

    result = gtd.add_pubchem_data(df, cfg)

    assert result.loc[0, "pubchem_cid"] == "Not Found"


def test_add_pubchem_data_reuses_resolution_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "canonical_smiles": ["C", "C"],
        }
    )
    cfg = pl.PubChemCfg(delay=0)

    calls: list[tuple[str | None, ...]] = []

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] | None = None,
        resolution_key: tuple[str | None, ...] | None = None,
        **kwargs: object,
    ) -> pl.PubChemResolution:
        if resolution_cache is not None and resolution_key is not None:
            if resolution_key in resolution_cache:
                return resolution_cache[resolution_key]
            calls.append(resolution_key)
            result = pl.PubChemResolution(cid="123", source="canonical_smiles")
            resolution_cache[resolution_key] = result
            return result
        raise AssertionError("resolution_cache should be provided")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)
    monkeypatch.setattr(
        pl,
        "get_properties",
        lambda cid, cfg: pl.Properties(None, None, None, None, None, None),
    )

    result = gtd.add_pubchem_data(df, cfg)

    assert len(calls) == 1
    assert result["pubchem_cid"].tolist() == ["123", "123"]


def test_add_pubchem_data_prefers_local_smiles(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame(
        {
            "canonical_smiles": ["C"],
            "pubchem_cid": ["1"],
            "pubchem_iupac_name": ["methane"],
            "pubchem_molecular_formula": ["CH4"],
            "pubchem_isomeric_smiles": ["C"],
            "pubchem_canonical_smiles": ["C"],
            "pubchem_inchi": ["InChI=1S/CH4/h1H4"],
            "pubchem_inchikey": ["VNWKTOKETHGBQD-UHFFFAOYSA-N"],
        }
    )
    cfg = pl.PubChemCfg(delay=0, prefer_local_smiles=True)

    def fake_resolve(*_: object, **__: object) -> pl.PubChemResolution:
        return pl.PubChemResolution(cid="321", source="cache")

    def fail(*_: object, **__: object) -> None:
        raise AssertionError("get_properties should not be called when prefer_local_smiles is enabled")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)
    monkeypatch.setattr(pl, "get_properties", fail)

    result = gtd.add_pubchem_data(df, cfg)

    expected = df.copy()
    for column in gtd.PUBCHEM_COLUMNS:
        expected[column] = expected[column].astype("string")
    pd.testing.assert_frame_equal(result, expected)



def test_add_pubchem_data_preserves_existing_values(

    monkeypatch: pytest.MonkeyPatch,
) -> None:
    df = pd.DataFrame(
        {

            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "canonical_smiles": ["", "C", "CC"],
            "molecule_type": ["Polymer", "Mixture", "Small molecule"],
            "pubchem_cid": ["111", "222", ""],
            "pubchem_iupac_name": ["existing", "mixture", ""],
        }
    )
    cfg = pl.PubChemCfg(
        delay=0,
        allow_polymer=False,
        prefer_local_smiles=False,
        prefer_local_values=True,
    )

    calls: list[str] = []
    warnings: list[tuple[str, dict[str, object]]] = []

    def fake_resolve(
        row: pd.Series,
        cache: dict[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        parent_loader: Callable[[str], pd.Series | None] | None = None,
        resolution_cache: dict[tuple[str | None, ...], pl.PubChemResolution] | None = None,
    ) -> str:
        calls.append(row["molecule_chembl_id"])
        return "333"

    monkeypatch.setattr(gtd, "resolve_pubchem_cid", fake_resolve)
    monkeypatch.setattr(pl, "get_properties", lambda cid, cfg: None)
    monkeypatch.setattr(
        gtd.logger,
        "warning",
        lambda event, **kwargs: warnings.append((event, kwargs)),
    )

    result = gtd.add_pubchem_data(df, cfg)

    assert calls == ["CHEMBL3"]
    assert any(
        event == "pubchem_skip_polymers"
        and kwargs.get("count") == 2
        and kwargs.get("polymer_count") == 1
        and kwargs.get("mixture_count") == 1
        and kwargs.get("indexes") == [0, 1]
        for event, kwargs in warnings
    )

    assert result.loc[0, "pubchem_cid"] == "111"
    assert result.loc[0, "pubchem_iupac_name"] == "existing"
    assert result.loc[1, "pubchem_cid"] == "222"
    assert result.loc[1, "pubchem_iupac_name"] == "mixture"
    assert result.loc[2, "pubchem_cid"] == "333"
 


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

    cid = gtd.resolve_pubchem_cid(row, cache, cfg, resolution_cache={})

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

    cid = gtd.resolve_pubchem_cid(row, cache, cfg, resolution_cache={})

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

    cid = gtd.resolve_pubchem_cid(
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

    monkeypatch.setattr(gtd.logger, "info", capture)
    monkeypatch.setattr(
        pl,
        "resolve_pubchem_record",
        lambda *args, **kwargs: pl.PubChemResolution(cid=None, source=None),
    )

    cid = gtd.resolve_pubchem_cid(
        row,
        cache,
        cfg,
        parent_loader=lambda _: None,
        resolution_cache={},
    )

    assert cid is None
    assert cache["CHEMBL2"] is None
    assert any(event == "pubchem_parent_structure_missing" for event, _ in events)


def test_add_pubchem_data_uses_disk_cache(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cache_path = tmp_path / "cid_cache.json"
    cache_path.write_text(json.dumps({"CHEMBL1": "321"}))

    df = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "canonical_smiles": ["C"],
        }
    )

    cfg = pl.PubChemCfg(delay=0, cid_cache_path=cache_path)

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg: pl.PubChemCfg,
        *,
        cid_cache: Mapping[str, str | None] | None = None,
        cache_key: str | None = None,
        **__: object,
    ) -> pl.PubChemResolution:
        assert cache_key == "CHEMBL1"
        assert cid_cache is not None and cid_cache.get("CHEMBL1") == "321"
        return pl.PubChemResolution(cid="321", source="cache")

    monkeypatch.setattr(pl, "resolve_pubchem_record", fake_resolve)

    props = pl.Properties("name", "formula", "i", "c", "inchi", "inchikey")
    monkeypatch.setattr(pl, "get_properties", lambda cid, cfg: props)

    result = gtd.add_pubchem_data(df, cfg)

    assert result.loc[0, "pubchem_cid"] == "321"
    assert result.loc[0, "pubchem_iupac_name"] == "name"
    cache_data = json.loads(cache_path.read_text())
    if "values" in cache_data:
        assert cache_data["values"] == {"CHEMBL1": "321"}
        assert cache_data["metadata"]["version"] == gtd._PUBCHEM_CACHE_SCHEMA_VERSION
    else:
        assert cache_data == {"CHEMBL1": "321"}


def test_write_pubchem_cid_cache_creates_parent_dir(tmp_path: Path) -> None:
    cache_path = tmp_path / "nested" / "cid_cache.json"

    assert not cache_path.parent.exists()

    gtd._write_pubchem_cid_cache(cache_path, {"CHEMBL1": "321"})

    assert cache_path.exists()
    payload = json.loads(cache_path.read_text())
    assert payload["values"] == {"CHEMBL1": "321"}
    assert payload["metadata"]["version"] == gtd._PUBCHEM_CACHE_SCHEMA_VERSION


def test_load_pubchem_cid_cache_expires(tmp_path: Path) -> None:
    cache_path = tmp_path / "cid_cache.json"
    expired_at = datetime.now(timezone.utc) - timedelta(hours=2)
    payload = {
        "metadata": {
            "version": gtd._PUBCHEM_CACHE_SCHEMA_VERSION,
            "updated_at": expired_at.isoformat(),
        },
        "values": {"CHEMBL1": "321"},
    }
    cache_path.write_text(json.dumps(payload))

    cache = gtd._load_pubchem_cid_cache(cache_path, ttl_hours=1)

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
    monkeypatch.setattr(
        gtd, "load_parent_catalog", lambda **__: {"CHEMBL1": "CHEMBL1_PARENT"}
    )
    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    monkeypatch.setattr(
        gtd,
        "add_pubchem_data",
        lambda df, cfg, **__: df,
    )
    captured_precomputed: dict[str, object] = {"data": None, "frame": None}

    def fake_attach_parent_molecule_ids(frame: pd.DataFrame, **kwargs: object):
        captured_precomputed["data"] = kwargs.get("precomputed")
        captured_precomputed["frame"] = frame.copy()
        return (

            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        )

    monkeypatch.setattr(gtd, "attach_parent_molecule_ids", fake_attach_parent_molecule_ids)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda p: "deadbeef")

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

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0

    prepared = captured_precomputed["data"]
    assert isinstance(prepared, gtd.ParentLookupPreparedData)
    captured_frame = captured_precomputed["frame"]
    assert isinstance(captured_frame, pd.DataFrame)
    expected_prepared = prepare_parent_lookup_data(
        captured_frame, cfg.molecule_catalog
    )
    pd.testing.assert_series_equal(
        prepared.child_ids, expected_prepared.child_ids, check_names=False
    )
    pd.testing.assert_series_equal(
        prepared.existing_parent_ids,
        expected_prepared.existing_parent_ids,
        check_names=False,
    )
    assert prepared.existing_parent_ids.tolist() == expected_prepared.existing_parent_ids.tolist()
    expected_need_lookup = set(
        expected_prepared.child_ids[
            (expected_prepared.child_ids != "")
            & (expected_prepared.existing_parent_ids == "")
        ]
    )
    assert prepared.need_lookup == expected_need_lookup

    available = set(captured.get("columns", []))
    expected_head = [c for c in TestitemsSchema.columns if c in available]
    expected_tail = sorted(available - set(TestitemsSchema.columns))
    assert captured["col_order"] == expected_head + expected_tail


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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, pubchem_cfg, **__: frame)
    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
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

    monkeypatch.setattr(gtd.pl, "init_session", fake_init_session)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda df, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert captured["init"] == (cfg.api, cfg.retry)


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
    monkeypatch.setattr(gtd.pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, *args, **__: frame)
    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(gtd, "query_parent_catalog", lambda *_, **__: {})
    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )
    monkeypatch.setattr(gtd, "normalize_testitems", lambda frame: frame)
    monkeypatch.setattr(gtd, "add_pipeline_metadata", lambda frame: frame)
    monkeypatch.setattr(
        gtd,
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
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "update_parent_catalog_cache", lambda *_, **__: None)

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert ids_source.iterations == 1
    assert received_ids is not None
    assert not isinstance(received_ids, list)
    assert received_ids.__class__.__name__ == "_tee"


def test_run_chembl_calls_pubchem_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "testitems.csv"
    input_csv.write_text(
        "molecule_chembl_id\nCHEMBL1\n", encoding=cfg.io.csv_encoding
    )

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
    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
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
    monkeypatch.setattr(gtd.pl, "init_session", lambda api, retry: None)
    monkeypatch.setattr(
        io,
        "write_csv",
        lambda frame, path, *, cfg, key_cols=None, col_order=None, **__: path,
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    rc = gtd.run_chembl(cfg, args)

    assert rc == 0
    assert len(resolve_calls) == 1


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

    monkeypatch.setattr(gtd, "query_parent_catalog", fake_query_parent_catalog)
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _, **__: frame)
    captured_catalog: dict[str, str] | None = None
    captured_source: str | None = None

    def fake_attach_parent_molecule_ids(
        frame: pd.DataFrame, **kwargs: object
    ) -> tuple[pd.DataFrame, gtd.ParentLookupStats]:
        nonlocal captured_catalog, captured_source
        captured_catalog = kwargs.get("catalog")
        captured_source = kwargs.get("source")

        catalog_cfg = kwargs.get("catalog_cfg")
        mapping = captured_catalog or {}
        parent_field = (
            catalog_cfg.parent_field if catalog_cfg is not None else "parent_molecule_chembl_id"
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

            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        )

    monkeypatch.setattr(gtd, "attach_parent_molecule_ids", fake_attach_parent_molecule_ids)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

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

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0
    assert captured_df is not None
    assert captured_df["parent_molecule_chembl_id"].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_EXISTING",
    ]
    assert query_calls == 1
    assert captured_catalog is parent_catalog
    assert captured_source == gtd.PARENT_LOOKUP_SOURCE_CACHE


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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _, **__: frame)

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
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    update_calls: list[dict[str, str]] = []
    original_update = gtd.update_parent_catalog_cache

    def tracking_update(catalog: Mapping[str, str], catalog_cfg: object) -> None:
        update_calls.append(dict(catalog))
        original_update(catalog, catalog_cfg)

    monkeypatch.setattr(gtd, "update_parent_catalog_cache", tracking_update)

    attach_calls: list[tuple[dict[str, str], str | None]] = []

    def fake_attach(
        frame: pd.DataFrame, **kwargs: object
    ) -> tuple[pd.DataFrame, gtd.ParentLookupStats]:
        catalog = dict(kwargs.get("catalog", {}))
        source_name = kwargs.get("source")
        attach_calls.append((catalog, source_name))
        stats_source = source_name or gtd.PARENT_LOOKUP_SOURCE_SKIPPED
        stats = gtd.ParentLookupStats(
            source=stats_source,
            missing=0,
            unique=len(catalog),
            attached=len(frame),
            uncovered=0,
        )
        return frame, stats

    monkeypatch.setattr(gtd, "attach_parent_molecule_ids", fake_attach)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

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
        gtd.PARENT_LOOKUP_SOURCE_PARTIAL,
        gtd.PARENT_LOOKUP_SOURCE_CACHE,
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

    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        lambda *_, **__: {},
    )

    result, stats = gtd.attach_parent_molecule_ids(
        source,
        client=object(),
        api_cfg=cfg.api,
        catalog_cfg=cfg.molecule_catalog,
        timeout=None,
    )

    expected = pd.Series(
        ["CHEMBL1_EXISTING", pd.NA], index=result.index, dtype="string"
    )
    pd.testing.assert_series_equal(
        result[parent_field], expected, check_names=False
    )
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
    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            ),
        ),
    )
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

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

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        gtd,
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

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _, **__: frame)
    monkeypatch.setattr(gtd, "analyze_table_quality", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "write_meta_yaml", lambda **kwargs: None)
    monkeypatch.setattr(gtd, "file_sha256", lambda path: "deadbeef")

    monkeypatch.setattr(
        gtd,
        "query_parent_catalog",
        lambda *_, **__: (_ for _ in ()).throw(requests.RequestException("boom")),
    )

    monkeypatch.setattr(
        gtd,
        "attach_parent_molecule_ids",
        lambda frame, **kwargs: (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
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

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    errors: list[tuple[str, dict[str, object]]] = []

    def fake_logger_error(event: str, *args: object, **kwargs: object) -> None:
        errors.append((event, kwargs))

    monkeypatch.setattr(gtd.logger, "error", fake_logger_error)

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
        gtd,
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
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
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
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_SYNC


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

    monkeypatch.setattr(gtd, "query_parent_catalog", fake_query)
    monkeypatch.setattr(gtd.molecule_catalog, "fetch_parent_catalog_for", fake_fetch)
    monkeypatch.setattr(gtd, "write_parent_catalog_cache", lambda *_, **__: None)
    monkeypatch.setattr(gtd, "update_parent_catalog_cache", lambda *_, **__: None)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
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
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_PARTIAL


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
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    monkeypatch.setattr(
        gtd,
        "load_parent_catalog",
        lambda **_: (_ for _ in ()).throw(
            AssertionError("load_parent_catalog should not be called")
        ),
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        precomputed=precomputed,
    )

    assert fetch_calls == [["CHEMBL10", "CHEMBL11"]]
    assert result[parent_field].tolist() == ["CHEMBL10_PARENT", "CHEMBL11_PARENT"]
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_PARTIAL


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

    base_catalog = {
        f"CHEMBL{i}": f"CHEMBL{i}_PARENT"
        for i in range(1, 5000)
    }

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
    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        catalog=tracking_catalog,
        source=gtd.PARENT_LOOKUP_SOURCE_CACHE,
        precomputed=precomputed,
    )

    expected_values = [f"CHEMBL{i}_PARENT" for i in range(1, 4)]
    assert result[catalog_cfg.parent_field].tolist() == expected_values
    assert stats.unique == 3
    assert stats.attached == 3
    assert stats.missing == 0
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_CACHE
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

    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})

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

    monkeypatch.setattr(gtd.molecule_catalog, "fetch_parent_catalog_for", fake_fetch)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
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
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_SYNC


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
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    first_result, first_stats = gtd.attach_parent_molecule_ids(
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
    assert first_stats.source == gtd.PARENT_LOOKUP_SOURCE_PARTIAL

    stored_catalog = json.loads(catalog_cfg.cache_path.read_text(encoding="utf-8"))
    assert stored_catalog == {
        "CHEMBL1": "CHEMBL1_PARENT",
        "CHEMBL2": "CHEMBL2_PARENT",
    }

    precomputed_second = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    second_result, second_stats = gtd.attach_parent_molecule_ids(
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
    assert second_stats.source == gtd.PARENT_LOOKUP_SOURCE_CACHE


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

    original_loads = gtd.molecule_catalog.json.loads
    calls = {"loads": 0}

    def counting_loads(data: str, *args: object, **kwargs: object) -> object:
        calls["loads"] += 1
        return original_loads(data, *args, **kwargs)

    monkeypatch.setattr(gtd.molecule_catalog.json, "loads", counting_loads)

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    first_result, _ = gtd.attach_parent_molecule_ids(
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
    second_result, _ = gtd.attach_parent_molecule_ids(
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
        gtd,
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
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        failing_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
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
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_SYNC


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

    monkeypatch.setattr(gtd, "load_parent_catalog", unexpected_load_parent_catalog)


    def unexpected_fetch(
        ids: list[str],
        *,
        client: object,
        api_cfg: object,
        timeout: float | None,
    ) -> dict[str, str]:
        raise AssertionError("fetch_parent_catalog_for should not be called")

    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        unexpected_fetch,
    )

    precomputed = (
        prepare_parent_lookup_data(df, catalog_cfg) if use_precomputed else None
    )
    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        catalog={"CHEMBL1": "CHEMBL1_PARENT"},
        source=gtd.PARENT_LOOKUP_SOURCE_CACHE,
        precomputed=precomputed,
    )

    assert result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]
    assert stats.missing == 0
    assert stats.attached == 1
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_CACHE

