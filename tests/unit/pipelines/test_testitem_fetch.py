from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import pandas as pd
import pytest

from library.config import ApiCfg, TestitemBatchRetryCfg
from library.pipelines.testitem import cli


class _FakeChemblLib:
    def __init__(self, frame: pd.DataFrame) -> None:
        self._frame = frame
        self.calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []

    def get_testitem(self, *args: Any, **kwargs: Any) -> pd.DataFrame:
        self.calls.append((args, kwargs))
        return self._frame.copy()


class _SentinelClient:
    """Client stub passed through to the ChemBL library."""

    def __repr__(self) -> str:  # pragma: no cover - debug helper
        return "<SentinelClient>"


def _collect(iterator: Iterator[pd.DataFrame]) -> list[pd.DataFrame]:
    return list(iterator)


@pytest.mark.unit
def test_fetch_testitems__does_not_reintroduce_nested_columns(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1"], dtype="string"),
            "pubchem_cid": pd.Series(["CID1"], dtype="string"),
            "standard_inchi": pd.Series(["InChI=1"], dtype="string"),
        }
    )
    fake_library = _FakeChemblLib(frame)
    monkeypatch.setattr(cli, "_load_chembl_library", lambda: fake_library)

    status, chunk_iter, requested = cli.fetch_testitems(
        ["CHEMBL1"],
        api_cfg=ApiCfg(),
        batch_size=5,
        timeout=1.0,
        client=_SentinelClient(),
        sample_ids=("CHEMBL1",),
        fields=(
            "molecule_structures.standard_inchi",
            "pubchem.cid",
        ),
        page_limit=10,
        retry_cfg=TestitemBatchRetryCfg(),
    )

    assert status == 0
    assert requested == ("CHEMBL1",)
    assert chunk_iter is not None

    chunks = _collect(chunk_iter)
    assert len(chunks) == 1
    chunk = chunks[0]

    assert set(chunk.columns) == {
        "molecule_chembl_id",
        "pubchem_cid",
        "standard_inchi",
    }
    assert chunk["pubchem_cid"].tolist() == ["CID1"]
    assert chunk["standard_inchi"].tolist() == ["InChI=1"]


@pytest.mark.unit
@pytest.mark.pipeline_scenario(
    "csv_loading",
    "normalization",
    "enrichment",
    "transformation_rules",
    "missing_data",
    "logging",
    "assembly",
    "export",
    "degradation",
    "idempotence",
)
def test_fetch_testitems__renames_pubchem_and_structure_columns(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL1"], dtype="string"),
            "molecule_structures.canonical_smiles": pd.Series(["C"], dtype="string"),
            "molecule_structures.standard_inchi": pd.Series(
                ["InChI=1S/C"], dtype="string"
            ),
            "molecule_structures.standard_inchi_key": pd.Series(
                ["XLYOFNOQVPJJNP-UHFFFAOYSA-N"], dtype="string"
            ),
            "pubchem.cid": pd.Series(["2244"], dtype="string"),
            "pubchem.canonical_smiles": pd.Series(["C"], dtype="string"),
            "pubchem.isomeric_smiles": pd.Series(["C"], dtype="string"),
            "pubchem.inchi": pd.Series(["InChI=1S/C"], dtype="string"),
            "pubchem.inchi_key": pd.Series(
                ["XLYOFNOQVPJJNP-UHFFFAOYSA-N"], dtype="string"
            ),
            "pubchem.iupac_name": pd.Series(["methane"], dtype="string"),
            "pubchem.molecular_formula": pd.Series(["CH4"], dtype="string"),
        }
    )
    fake_library = _FakeChemblLib(frame)
    monkeypatch.setattr(cli, "_load_chembl_library", lambda: fake_library)

    fields = (
        "molecule_structures.canonical_smiles",
        "molecule_structures.standard_inchi",
        "molecule_structures.standard_inchi_key",
        "pubchem.cid",
        "pubchem.canonical_smiles",
        "pubchem.isomeric_smiles",
        "pubchem.inchi",
        "pubchem.inchi_key",
        "pubchem.iupac_name",
        "pubchem.molecular_formula",
    )

    status, chunk_iter, _ = cli.fetch_testitems(
        ["CHEMBL1"],
        api_cfg=ApiCfg(),
        batch_size=5,
        timeout=1.0,
        client=_SentinelClient(),
        sample_ids=("CHEMBL1",),
        fields=fields,
        page_limit=10,
        retry_cfg=TestitemBatchRetryCfg(),
    )

    assert status == 0
    assert chunk_iter is not None

    chunks = _collect(chunk_iter)
    assert len(chunks) == 1
    chunk = chunks[0]

    expected_columns = {
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_inchi",
        "standard_inchi_key",
        "pubchem_cid",
        "pubchem_canonical_smiles",
        "pubchem_isomeric_smiles",
        "pubchem_inchi",
        "pubchem_inchikey",
        "pubchem_iupac_name",
        "pubchem_molecular_formula",
    }

    assert expected_columns.issubset(set(chunk.columns))
    assert chunk.loc[0, "canonical_smiles"] == "C"
    assert chunk.loc[0, "standard_inchi"] == "InChI=1S/C"
    assert chunk.loc[0, "standard_inchi_key"] == "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    assert chunk.loc[0, "pubchem_cid"] == "2244"
    assert chunk.loc[0, "pubchem_canonical_smiles"] == "C"
    assert chunk.loc[0, "pubchem_inchi"] == "InChI=1S/C"
    assert chunk.loc[0, "pubchem_inchikey"] == "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
    assert chunk.loc[0, "pubchem_iupac_name"] == "methane"
    assert chunk.loc[0, "pubchem_molecular_formula"] == "CH4"
