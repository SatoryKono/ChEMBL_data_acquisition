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

