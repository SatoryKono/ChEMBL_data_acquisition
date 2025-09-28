from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library.processing import chembl_library as cl
from library.clients import pubchem_library as pl
from library import testitem_enrichment
from scripts import get_testitem_data


def _cleanup_output(path: Path) -> None:
    for candidate in (
        path,
        path.with_name(path.name + ".meta.yaml"),
        path.with_suffix(".quality.json"),
        path.with_name(f"{path.stem}_failure_cases.csv"),
    ):
        if candidate.exists() and candidate.is_file():
            candidate.unlink()


@pytest.mark.usefixtures("smoke_output_dir")
def test_testitem_salt_flags_smoke(
    smoke_output_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = Path("tests/data/input-smoke/testitem.csv")
    output_csv = smoke_output_dir / "testitem_flags.csv"
    _cleanup_output(output_csv)

    def fake_get_testitem(ids, cfg, client, chunk_size, timeout):  # type: ignore[no-untyped-def]
        rows: list[dict[str, object]] = []
        for mol_id in ids:
            if mol_id == "CHEMBL1":
                child = "CHEMBL100"
                parent = None
            elif mol_id == "CHEMBL2":
                child = "CHEMBL200"
                parent = "CHEMBL0200P"
            else:
                child = str(mol_id)
                parent = None
            rows.append(
                {
                    "molecule_chembl_id": child,
                    "parent_molecule_chembl_id": parent,
                    "max_phase": "4",
                }
            )
        rows.append(
            {
                "molecule_chembl_id": "CHEMBL300",
                "parent_molecule_chembl_id": "CHEMBL300",
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

    hierarchy = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL100", "CHEMBL200", "CHEMBL201"],
            "parent_molecule_chembl_id": ["CHEMBL0100P", "CHEMBL0200P", "CHEMBL0201P"],
        }
    )
    catalog = pd.DataFrame(
        {
            "molecule_chembl_id": [
                "CHEMBL100",
                "CHEMBL0100P",
                "CHEMBL200",
                "CHEMBL0200P",
                "CHEMBL300",
            ],
            "natural_product": ["Y", "N", "N", "", ""],
            "prodrug": ["", "Y", "N", "N", ""],
            "polymer_flag": ["0", "0", "1", "0", ""],
        }
    )

    monkeypatch.setattr(cl, "get_testitem", fake_get_testitem)
    monkeypatch.setattr(pl, "init_session", lambda *_, **__: None)
    monkeypatch.setattr(pl, "get_cid_from_smiles", fake_get_cid)
    monkeypatch.setattr(pl, "get_properties", fake_get_properties)
    monkeypatch.setattr(
        testitem_enrichment,
        "_load_sources",
        lambda cfg, io_cfg: testitem_enrichment._SourceFrames(
            hierarchy=hierarchy, catalog=catalog
        ),
    )
    monkeypatch.setattr(
        get_testitem_data, "analyze_table_quality", lambda *_, **__: None
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
            str(Path("config.yaml")),
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()

    df = pd.read_csv(output_csv)
    assert {"salt_chembl_id", "natural_product"}.issubset(df.columns)

    salt_row = df.loc[df["molecule_chembl_id"] == "CHEMBL100"].iloc[0]
    assert salt_row["salt_chembl_id"] == "CHEMBL100"
    assert salt_row["prodrug"] is True

    direct_row = df.loc[df["molecule_chembl_id"] == "CHEMBL200"].iloc[0]
    assert direct_row["salt_chembl_id"] == "CHEMBL200"
    assert direct_row["natural_product"] is False

    nonsalt_row = df.loc[df["molecule_chembl_id"] == "CHEMBL300"].iloc[0]
    assert pd.isna(nonsalt_row["salt_chembl_id"])
