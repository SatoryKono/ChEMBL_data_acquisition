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

    class FakeCatalog:
        def __init__(self, mapping: dict[str, str]) -> None:
            self._mapping = mapping
            self.was_refreshed = False

        def __bool__(self) -> bool:
            return bool(self._mapping)

        def lookup(self, children: list[str]) -> dict[str, str]:
            return {
                child: self._mapping[child]
                for child in children
                if child in self._mapping
            }

    def fake_load_parent_catalog(**_: object) -> FakeCatalog:
        frame = pd.read_csv(catalog_csv)
        mapping = {
            str(row["molecule_chembl_id"]): str(row["parent_molecule_chembl_id"])
            for _, row in frame.iterrows()
        }
        return FakeCatalog(mapping)

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(get_testitem_data, "load_parent_catalog", fake_load_parent_catalog)
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
    assert list(df["parent_molecule_chembl_id"]) == [
        "CHEMBL9001",
        "CHEMBL9002",
    ]
