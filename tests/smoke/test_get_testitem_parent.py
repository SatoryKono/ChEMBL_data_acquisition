from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library import chembl_library as cl
from library import pubchem_library as pl
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


def test_get_testitem_parent_catalog(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Smoke-test parent catalogue enrichment for test items."""

    input_csv = Path("tests/data/input-smoke/testitem.csv")
    catalog_csv = Path("tests/data/input-smoke/molecule_catalog.csv")
    output_csv = smoke_output_dir / "testitem_parent.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
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

    def fake_load_parent_catalog(**_: object) -> dict[str, str]:
        frame = pd.read_csv(catalog_csv)
        return {
            str(row["molecule_chembl_id"]): str(row["parent_molecule_chembl_id"])
            for _, row in frame.iterrows()
        }

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(get_testitem_data, "load_parent_catalog", fake_load_parent_catalog)
    monkeypatch.setattr(
        get_testitem_data,
        "_cache_state",
        lambda _path: (True, 0.0),
    )
    fetch_calls: list[list[str]] = []

    def fake_fetch_parent_catalog_for(ids, **_: object) -> dict[str, str]:
        fetch_calls.append(list(ids))
        return {}

    monkeypatch.setattr(
        get_testitem_data.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setattr(get_testitem_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
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
    assert fetch_calls == []


def test_get_testitem_skips_parent_lookup_when_present(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure parent lookups are skipped when data already contains parents."""

    input_csv = Path("tests/data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem_parent_present.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
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

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(get_testitem_data, "load_parent_catalog", fail_load_parent_catalog)
    monkeypatch.setattr(
        get_testitem_data.molecule_catalog,
        "fetch_parent_catalog_for",
        fake_fetch_parent_catalog_for,
    )
    monkeypatch.setattr(get_testitem_data, "analyze_table_quality", lambda *_, **__: None)

    exit_code = get_testitem_data.main(
        [
            "--input",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
            "--config",
            str(Path("config.yaml")),
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert list(df["parent_molecule_chembl_id"]) == ["CHEMBL9001", "CHEMBL9002"]
    assert fetch_called is False
