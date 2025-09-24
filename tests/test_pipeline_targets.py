"""Tests for :mod:`library.pipeline_targets`."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Any

import pandas as pd

from library.pipeline_targets import PipelineResult, run_pipeline


def _iterator() -> Iterator[Iterable[str]]:
    yield ["CHEMBL1", "CHEMBL2"]


def test_run_pipeline_passes_batch_size() -> None:
    """Ensure ``batch_size`` is forwarded to fetchers that accept it."""

    received: dict[str, Any] = {}

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]],
        cfg: object,
        *,
        batch_size: int,
    ) -> pd.DataFrame:
        received["batch_size"] = batch_size
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
        batch_size=25,
    )

    assert isinstance(result, PipelineResult)
    assert received["batch_size"] == 25
    assert result.chembl.loc[0, "target_chembl_id"] == "CHEMBL1"


def test_run_pipeline_ignores_optional_keywords_for_wrappers() -> None:
    """Wrappers without ``batch_size`` support are retried without the keyword."""

    calls = {"count": 0}

    def base_fetcher(chunks: Iterator[Iterable[str]], cfg: object) -> pd.DataFrame:
        calls["count"] += 1
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    def cached_fetcher(chunks: Iterator[Iterable[str]], cfg: object) -> pd.DataFrame:
        # Simulate a closure that intentionally omits ``batch_size``.
        return base_fetcher(chunks, cfg)

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=cached_fetcher,
        batch_size=15,
    )

    assert calls["count"] == 1
    assert list(result.chembl["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]


def test_optional_stages_receive_dataframe() -> None:
    """Optional stages operate on the ChEMBL output."""

    def chembl_fetcher(
        chunks: Iterator[Iterable[str]],
        cfg: object,
        *,
        batch_size: int,
    ) -> pd.DataFrame:
        ids = [item for chunk in chunks for item in chunk]
        return pd.DataFrame({"target_chembl_id": ids})

    def uniprot_fetcher(
        df: pd.DataFrame, cfg: object | None = None, **_: Any
    ) -> pd.DataFrame:
        assert list(df["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
        return pd.DataFrame({"uniprot_id": ["P1", "P2"]})

    result = run_pipeline(
        _iterator,
        chembl_cfg=object(),
        chembl_fetcher=chembl_fetcher,
        uniprot_fetcher=uniprot_fetcher,
        uniprot_cfg=object(),
    )

    assert result.uniprot is not None
    assert list(result.uniprot["uniprot_id"]) == ["P1", "P2"]

    # Mapping interface exposes the stored frames.
    assert result["chembl"].equals(result.chembl)
    assert result.as_dict()["uniprot"].equals(result.uniprot)
