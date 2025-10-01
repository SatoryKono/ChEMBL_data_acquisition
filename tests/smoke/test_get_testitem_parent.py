from __future__ import annotations

import csv
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library import chembl_library as cl
from library import cli as base_cli
from library import pubchem_library as pl
from library.clients import pubchem as pc
from library.config import IoCfg
from library.utils.config import DEFAULT_CONFIG_RELATIVE

_ORIGINAL_APPLY = base_cli.apply_config_overrides
from scripts import get_testitem_data


def _cleanup_output(path: Path) -> None:
    """Remove existing artefacts to keep assertions deterministic."""

    for candidate in (
        path,
        path.with_name(path.name + ".meta.yaml"),
        path.with_suffix(".quality.json"),
        path.with_name(f"{path.stem}_failure_cases.csv"),
    ):
        if candidate.exists() and candidate.is_file():
            candidate.unlink()


def test_read_input_ids_handles_large_input(tmp_path: Path) -> None:
    """Ensure identifiers are collected once even for large inputs."""

    total_ids = 1500
    input_csv = tmp_path / "testitem_large.csv"
    with input_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["chembl_id"])
        for index in range(total_ids):
            writer.writerow([f"CHEMBL{index}"])

    status, result = get_testitem_data.read_input_ids(
        input_csv,
        column="chembl_id",
        io_cfg=IoCfg(),
        limit=None,
    )

    assert status == 0
    assert result is not None
    assert isinstance(result.ids, list)
    assert len(result.ids) == total_ids
    assert result.requested_ids is result.ids
    assert list(result.ids_iter) == result.ids
    expected_sample_size = min(
        total_ids,
        get_testitem_data._FETCH_ERROR_SAMPLE_SIZE,
    )
    assert result.sample_ids == tuple(result.ids[:expected_sample_size])


def test_get_testitem_parent_catalog(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test parent catalogue enrichment for test items."""

    input_csv = Path("tests/data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem_parent.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, mol_id in enumerate(ids, start=1):
            rows.append(
                {
                    "molecule_chembl_id": str(mol_id),
                    "canonical_smiles": f"C{idx}H{2 * idx + 2}",
                    "max_phase": "4",
                }
            )
        return pd.DataFrame(rows)

    def fake_get_cid(smiles: str, cfg):  # type: ignore[no-untyped-def]
        return "12345"

    def fake_get_properties(cid: str, cfg):  # type: ignore[no-untyped-def]
        return SimpleNamespace(
            IUPACName="Mock acid",
            MolecularFormula="C2H6O",
            iSMILES="CCO",
            cSMILES="CCO",
            InChI="InChI=1S/C2H6O",
            InChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )

    def fail_load_parent_catalog(**_: object) -> dict[str, str]:  # pragma: no cover
        pytest.fail("load_parent_catalog should not be triggered for partial coverage")

    original_apply = _ORIGINAL_APPLY

    def patched_apply(*args, **kwargs):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.pubchem.resolve_order = ("cache", "smiles")
        return cfg

    monkeypatch.setattr(
        get_testitem_data.cli, "apply_config_overrides", patched_apply
    )
    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        get_testitem_data, "load_parent_catalog", fail_load_parent_catalog
    )
    def fake_query_parent_catalog(*args: object, **kwargs: object) -> dict[str, str]:  # type: ignore[no-untyped-def]
        return {}

    monkeypatch.setattr(
        get_testitem_data,
        "query_parent_catalog",
        fake_query_parent_catalog,
    )
    fetch_calls: list[list[str]] = []
    mapping = {
        "CHEMBL1": "CHEMBL9001",
        "CHEMBL2": "CHEMBL9002",
    }

    def fake_fetch_parent_catalog_for(ids, **_: object) -> dict[str, str]:
        fetch_calls.append(list(ids))
        return {key: mapping[key] for key in ids if key in mapping}

    monkeypatch.setattr(
        get_testitem_data.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    original_apply = _ORIGINAL_APPLY

    def patched_apply_config_overrides(*args: object, **kwargs: object):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.pubchem.resolve_order = ("smiles",)
        return cfg

    monkeypatch.setattr(
        get_testitem_data,
        "write_parent_catalog_cache",
        lambda data, cfg: None,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "update_parent_catalog_cache",
        lambda data, cfg: None,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "load_molecule_hierarchy_lookup",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        get_testitem_data.prepare_parent_enrichment.__globals__["molecule_catalog"],
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setitem(
        get_testitem_data.prepare_parent_enrichment.__globals__,
        "query_parent_catalog",
        fake_query_parent_catalog,
    )
    monkeypatch.setattr(
        get_testitem_data.run_parent_enrichment.__globals__["molecule_catalog"],
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setitem(
        get_testitem_data.run_parent_enrichment.__globals__,
        "query_parent_catalog",
        fake_query_parent_catalog,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "apply_testitem_enrichment",
        lambda frame, *, enrichment_cfg, io_cfg: (0, frame),
    )
    def fake_run_parent_enrichment(prep, **kwargs):  # type: ignore[no-untyped-def]
        df = prep.df.copy()
        df["parent_molecule_chembl_id"] = (
            df["molecule_chembl_id"].map(mapping).astype("string")
        )
        stats = get_testitem_data.ParentLookupStats(
            source=get_testitem_data.PARENT_LOOKUP_SOURCE_LOOKUP,
            missing=0,
            unique=len(df["molecule_chembl_id"].unique()),
            attached=len(df),
            uncovered=0,
        )
        return 0, get_testitem_data.ParentEnrichmentResult(df=df, parent_stats=stats)

    monkeypatch.setattr(
        get_testitem_data,
        "run_parent_enrichment",
        fake_run_parent_enrichment,
    )
    monkeypatch.setattr(
        get_testitem_data, "analyze_table_quality", lambda *_, **__: None
    )
    monkeypatch.setattr(
        get_testitem_data.cli,
        "apply_config_overrides",
        patched_apply_config_overrides,
    )

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            CONFIG_CLI_PATH,
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert "parent_molecule_chembl_id" in df.columns
    assert df["parent_molecule_chembl_id"].isna().sum() == 0
    assert list(df["parent_molecule_chembl_id"]) == [
        "CHEMBL9001",
        "CHEMBL9002",
    ]
    assert fetch_calls and len(fetch_calls) == 1
    assert set(fetch_calls[0]) == {"CHEMBL1", "CHEMBL2"}


def test_get_testitem_skips_parent_lookup_when_present(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure parent lookups are skipped when data already contains parents."""

    input_csv = Path("tests/data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem_parent_present.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, mol_id in enumerate(ids, start=1):
            rows.append(
                {
                    "molecule_chembl_id": str(mol_id),
                    "parent_molecule_chembl_id": f"CHEMBL90{idx:02d}",
                    "canonical_smiles": f"C{idx}H{2 * idx + 2}",
                    "max_phase": "4",
                }
            )
        return pd.DataFrame(rows)

    def fake_get_cid(smiles: str, cfg):  # type: ignore[no-untyped-def]
        return "12345"

    def fake_get_properties(cid: str, cfg):  # type: ignore[no-untyped-def]
        return SimpleNamespace(
            IUPACName="Mock acid",
            MolecularFormula="C2H6O",
            iSMILES="CCO",
            cSMILES="CCO",
            InChI="InChI=1S/C2H6O",
            InChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )

    def fail_load_parent_catalog(**_: object) -> dict[str, str]:
        pytest.fail("load_parent_catalog should not be called when parents exist")

    fetch_called = False

    def fake_fetch_parent_catalog_for(ids, **_: object) -> dict[str, str]:
        nonlocal fetch_called
        fetch_called = True
        return {}

    original_apply = _ORIGINAL_APPLY

    def patched_apply(*args, **kwargs):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.pubchem.resolve_order = ("cache", "smiles")
        return cfg

    monkeypatch.setattr(
        get_testitem_data.cli, "apply_config_overrides", patched_apply
    )
    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        get_testitem_data, "load_parent_catalog", fail_load_parent_catalog
    )
    monkeypatch.setattr(
        get_testitem_data.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    original_apply = _ORIGINAL_APPLY

    def patched_apply_config_overrides(*args: object, **kwargs: object):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.pubchem.resolve_order = ("smiles",)
        return cfg

    monkeypatch.setattr(
        get_testitem_data, "analyze_table_quality", lambda *_, **__: None
    )
    monkeypatch.setattr(
        get_testitem_data.cli,
        "apply_config_overrides",
        patched_apply_config_overrides,
    )

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            CONFIG_CLI_PATH,
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert list(df["parent_molecule_chembl_id"]) == ["CHEMBL9001", "CHEMBL9002"]
    assert fetch_called is False


def test_get_testitem_refreshes_outdated_parents(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure existing parent IDs are refreshed when forced via configuration."""

    input_csv = Path("tests/data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem_parent_refresh.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout, **kwargs):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for idx, mol_id in enumerate(ids, start=1):
            rows.append(
                {
                    "molecule_chembl_id": str(mol_id),
                    "parent_molecule_chembl_id": f"CHEMBL99{idx:02d}",
                    "canonical_smiles": f"C{idx}H{2 * idx + 2}",
                    "max_phase": "4",
                }
            )
        return pd.DataFrame(rows)

    def fake_get_cid(smiles: str, cfg):  # type: ignore[no-untyped-def]
        return "12345"

    def fake_get_properties(cid: str, cfg):  # type: ignore[no-untyped-def]
        return SimpleNamespace(
            IUPACName="Mock acid",
            MolecularFormula="C2H6O",
            iSMILES="CCO",
            cSMILES="CCO",
            InChI="InChI=1S/C2H6O",
            InChIKey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        )

    mapping = {
        "CHEMBL1": "CHEMBL9001",
        "CHEMBL2": "CHEMBL9002",
    }
    fetch_calls: list[list[str]] = []

    def fake_fetch_parent_catalog_for(ids, **_: object) -> dict[str, str]:
        fetch_calls.append(list(ids))
        return {key: mapping[key] for key in ids if key in mapping}

    original_apply = _ORIGINAL_APPLY

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pc, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        get_testitem_data.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "write_parent_catalog_cache",
        lambda data, cfg: None,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "update_parent_catalog_cache",
        lambda data, cfg: None,
    )
    monkeypatch.setattr(
        get_testitem_data,
        "load_molecule_hierarchy_lookup",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        get_testitem_data, "analyze_table_quality", lambda *_, **__: None
    )
    monkeypatch.setattr(get_testitem_data, "query_parent_catalog", lambda *_, **__: {})

    def patched_apply_config_overrides(*args: object, **kwargs: object):  # type: ignore[no-untyped-def]
        cfg = original_apply(*args, **kwargs)
        cfg.sources.chembl.molecule_catalog.force_refresh_existing = True
        cfg.pubchem.resolve_order = ("smiles",)
        return cfg

    monkeypatch.setattr(
        get_testitem_data.cli,
        "apply_config_overrides",
        patched_apply_config_overrides,
    )

    monkeypatch.setattr(
        get_testitem_data.prepare_parent_enrichment.__globals__["molecule_catalog"],
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setitem(
        get_testitem_data.prepare_parent_enrichment.__globals__,
        "query_parent_catalog",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        get_testitem_data.run_parent_enrichment.__globals__["molecule_catalog"],
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setitem(
        get_testitem_data.run_parent_enrichment.__globals__,
        "query_parent_catalog",
        lambda *_, **__: {},
    )
    monkeypatch.setattr(
        get_testitem_data,
        "apply_testitem_enrichment",
        lambda frame, *, enrichment_cfg, io_cfg: (0, frame),
    )

    def fake_run_parent_enrichment(prep, **kwargs):  # type: ignore[no-untyped-def]
        df = prep.df.copy()
        df["parent_molecule_chembl_id"] = (
            df["molecule_chembl_id"].map(mapping).astype("string")
        )
        stats = get_testitem_data.ParentLookupStats(
            source=get_testitem_data.PARENT_LOOKUP_SOURCE_LOOKUP,
            missing=0,
            unique=len(df["molecule_chembl_id"].unique()),
            attached=len(df),
            uncovered=0,
        )
        return 0, get_testitem_data.ParentEnrichmentResult(df=df, parent_stats=stats)

    monkeypatch.setattr(
        get_testitem_data,
        "run_parent_enrichment",
        fake_run_parent_enrichment,
    )

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            CONFIG_CLI_PATH,
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert list(df["parent_molecule_chembl_id"]) == ["CHEMBL9001", "CHEMBL9002"]
    assert fetch_calls and len(fetch_calls) == 1
    assert set(fetch_calls[0]) == {"CHEMBL1", "CHEMBL2"}
CONFIG_CLI_PATH = str(DEFAULT_CONFIG_RELATIVE)

