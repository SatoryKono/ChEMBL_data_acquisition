from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import pytest
import requests

from library import chembl_library as cl
from library import io
from library import pubchem_library as pl
from library.config import Config
from schemas import TestitemsSchema
from scripts import get_testitem_data as gtd


def test_add_pubchem_data_missing_uses_na(monkeypatch: pytest.MonkeyPatch) -> None:
    df = pd.DataFrame({"canonical_smiles": ["C"]})
    cfg = pl.PubChemCfg(delay=0)

    monkeypatch.setattr(pl, "get_cid_from_smiles", lambda *_: None)

    result = gtd.add_pubchem_data(df, cfg)
    pubchem_cols = [col for col in result.columns if col.startswith("pubchem_")]
    assert pubchem_cols
    assert result[pubchem_cols].isna().all().all()
    assert not (result[pubchem_cols] == "Not Found").any().any()


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

    def fail(*_: object, **__: object) -> None:  # pragma: no cover - defensive
        raise AssertionError("PubChem lookup should not be called")

    monkeypatch.setattr(pl, "get_cid_from_smiles", fail)
    monkeypatch.setattr(pl, "get_properties", fail)

    result = gtd.add_pubchem_data(df, cfg)

    expected = df.copy()
    for column in gtd.PUBCHEM_COLUMNS:
        expected[column] = expected[column].astype("string")
    pd.testing.assert_frame_equal(result, expected)


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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda df, cfg: df)
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
            ),
        ),
    )
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
        return output

    monkeypatch.setattr(io, "write_csv", fake_write_csv)

    rc = gtd.run_chembl(cfg, args)
    assert rc == 0

    available = set(df.columns) | {"pipeline_version", "timestamp_utc"}
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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, pubchem_cfg: frame)
    monkeypatch.setattr(gtd, "load_parent_catalog", lambda **__: {})

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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
    captured_catalog: dict[str, str] | None = None
    captured_source: str | None = None

    def fake_attach_parent_molecule_ids(
        frame: pd.DataFrame, **kwargs: object
    ) -> tuple[pd.DataFrame, gtd.ParentLookupStats]:
        nonlocal captured_catalog, captured_source
        captured_catalog = kwargs.get("catalog")
        captured_source = kwargs.get("source")
        return (
            frame,
            gtd.ParentLookupStats(
                source=gtd.PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
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
    monkeypatch.setattr(gtd, "add_pubchem_data", lambda frame, _: frame)
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


def test_attach_parent_molecule_ids_fetches_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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
    ) -> dict[str, str]:
        captured_ids["ids"] = list(ids)
        return fetched

    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    assert captured_ids["ids"] == ["CHEMBL2"]
    assert result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert stats.unique == 2
    assert stats.attached == 2
    assert stats.missing == 0
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_REMOTE


def test_attach_parent_molecule_ids_updates_cache_for_reuse(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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
    ) -> dict[str, str]:
        fetch_calls.append(list(ids))
        return {"CHEMBL2": "CHEMBL2_PARENT"}

    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch,
    )

    first_result, first_stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    assert fetch_calls == [["CHEMBL2"]]
    assert first_result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert first_stats.source == gtd.PARENT_LOOKUP_SOURCE_REMOTE

    stored_catalog = json.loads(catalog_cfg.cache_path.read_text(encoding="utf-8"))
    assert stored_catalog == {
        "CHEMBL1": "CHEMBL1_PARENT",
        "CHEMBL2": "CHEMBL2_PARENT",
    }

    second_result, second_stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    assert fetch_calls == [["CHEMBL2"]]
    assert second_result[catalog_cfg.parent_field].tolist() == [
        "CHEMBL1_PARENT",
        "CHEMBL2_PARENT",
    ]
    assert second_stats.source == gtd.PARENT_LOOKUP_SOURCE_CACHE


def test_attach_parent_molecule_ids_uses_sqlite_after_migration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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

    first_result, _ = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    assert calls["loads"] == 1
    assert first_result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]

    calls["loads"] = 0

    second_result, _ = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    assert calls["loads"] == 0
    assert second_result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]


def test_attach_parent_molecule_ids_fetch_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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
    ) -> dict[str, str]:
        raise requests.RequestException("boom")

    monkeypatch.setattr(
        gtd.molecule_catalog,
        "fetch_parent_catalog_for",
        failing_fetch,
    )

    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
    )

    parent_values = result[catalog_cfg.parent_field].tolist()
    assert parent_values == [pd.NA]
    assert stats.unique == 1
    assert stats.attached == 0
    assert stats.missing == 1
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_REMOTE


def test_attach_parent_molecule_ids_uses_cache_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
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

    result, stats = gtd.attach_parent_molecule_ids(
        df,
        client=object(),
        api_cfg=cfg.sources.chembl.api,
        catalog_cfg=catalog_cfg,
        timeout=None,
        catalog={"CHEMBL1": "CHEMBL1_PARENT"},
        source=gtd.PARENT_LOOKUP_SOURCE_CACHE,
    )

    assert result[catalog_cfg.parent_field].tolist() == ["CHEMBL1_PARENT"]
    assert stats.missing == 0
    assert stats.attached == 1
    assert stats.source == gtd.PARENT_LOOKUP_SOURCE_CACHE

