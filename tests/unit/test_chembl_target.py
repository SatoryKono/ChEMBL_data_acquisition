"""Unit tests for the ChEMBL target helpers."""

from __future__ import annotations

from collections.abc import Sequence

import pytest
import requests

from library.config import ApiCfg, UniprotMappingCfg
from library.pipelines.target.chembl_target import (
    TARGET_INCLUDE_PARAMS,
    iter_target_batches,
)


class _StubChemblClient:
    """Deterministic stub that mimics :class:`ChemblClient`."""

    def __init__(self, responses: dict[str, object]) -> None:
        self._responses = responses
        self.calls: list[tuple[str, float | None]] = []

    def request_json(
        self, url: str, *, cfg: ApiCfg, timeout: float | None = None
    ) -> dict:
        del cfg  # Unused in the stub.
        self.calls.append((url, timeout))
        response = self._responses[url]
        if isinstance(response, Exception):
            raise response
        assert isinstance(response, dict)
        return response


def _build_response(target_id: str, pref_name: str) -> dict[str, object]:
    return {
        "targets": [
            {
                "target_chembl_id": target_id,
                "pref_name": pref_name,
                "target_type": "SINGLE PROTEIN",
                "tax_id": 9606,
                "species_group_flag": 1,
                "target_components": [],
                "protein_classifications": [],
                "cross_references": [],
            }
        ]
    }


def _chunk_url(cfg: ApiCfg, ids: Sequence[str]) -> str:
    base = cfg.chembl_base.rstrip("/")
    joined_ids = ",".join(ids)
    return (
        f"{base}/target.json?format=json"
        f"&include={TARGET_INCLUDE_PARAMS}"
        f"&target_chembl_id__in={joined_ids}"
    )


@pytest.mark.unit
def test_iter_target_batches__splits_chunk_on_timeout(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """The iterator should split large chunks when a read timeout occurs."""

    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=12.0)
    mapping_cfg = UniprotMappingCfg()
    timeout = 9.5

    combined_url = _chunk_url(cfg, ["CHEMBL1", "CHEMBL2"])
    responses = {
        combined_url: requests.ReadTimeout("simulated timeout"),
        _chunk_url(cfg, ["CHEMBL1"]): _build_response("CHEMBL1", "Alpha"),
        _chunk_url(cfg, ["CHEMBL2"]): _build_response("CHEMBL2", "Beta"),
    }
    client = _StubChemblClient(responses)

    with caplog.at_level("WARNING"):
        batches = list(
            iter_target_batches(
                ["CHEMBL1", "CHEMBL2"],
                cfg=cfg,
                client=client,
                mapping_cfg=mapping_cfg,
                chunk_size=2,
                timeout=timeout,
            )
        )

    assert [call[0] for call in client.calls] == [
        combined_url,
        _chunk_url(cfg, ["CHEMBL1"]),
        _chunk_url(cfg, ["CHEMBL2"]),
    ]
    assert all(call[1] == timeout for call in client.calls)

    assert len(batches) == 2
    parsed_ids = [parsed[2]["target_chembl_id"].iat[0] for parsed in batches]
    assert parsed_ids == ["CHEMBL1", "CHEMBL2"]

    assert any(
        record.getMessage().startswith("chembl_timeout_split")
        for record in caplog.records
    )


@pytest.mark.unit
def test_iter_target_batches__propagates_timeout_without_split() -> None:
    """Disable fallback splitting so that higher level retry logic can react."""

    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=8.0)
    mapping_cfg = UniprotMappingCfg()
    timeout = 6.0

    combined_url = _chunk_url(cfg, ["CHEMBL10", "CHEMBL11"])
    responses = {combined_url: requests.ReadTimeout("simulated timeout")}
    client = _StubChemblClient(responses)

    with pytest.raises(requests.ReadTimeout):
        list(
            iter_target_batches(
                ["CHEMBL10", "CHEMBL11"],
                cfg=cfg,
                client=client,
                mapping_cfg=mapping_cfg,
                chunk_size=2,
                timeout=timeout,
                enable_split_fallback=False,
            )
        )

    assert client.calls == [(combined_url, timeout)]


@pytest.mark.unit
def test_iter_target_batches__splits_chunk_on_connection_error(
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Chunk-level connection errors should fall back to per-ID requests."""

    cfg = ApiCfg(chembl_base="https://example.test/api", timeout_read=7.0)
    mapping_cfg = UniprotMappingCfg()
    timeout = 5.0

    combined_url = _chunk_url(cfg, ["CHEMBL1", "CHEMBL2"])
    responses = {
        combined_url: requests.ConnectionError("simulated connection reset"),
        _chunk_url(cfg, ["CHEMBL1"]): _build_response("CHEMBL1", "Alpha"),
        _chunk_url(cfg, ["CHEMBL2"]): _build_response("CHEMBL2", "Beta"),
    }
    client = _StubChemblClient(responses)

    with caplog.at_level("WARNING"):
        batches = list(
            iter_target_batches(
                ["CHEMBL1", "CHEMBL2"],
                cfg=cfg,
                client=client,
                mapping_cfg=mapping_cfg,
                chunk_size=2,
                timeout=timeout,
            )
        )

    assert [call[0] for call in client.calls] == [
        combined_url,
        _chunk_url(cfg, ["CHEMBL1"]),
        _chunk_url(cfg, ["CHEMBL2"]),
    ]
    assert all(call[1] == timeout for call in client.calls)

    assert len(batches) == 2
    parsed_ids = [parsed[2]["target_chembl_id"].iat[0] for parsed in batches]
    assert parsed_ids == ["CHEMBL1", "CHEMBL2"]

    assert any(
        record.getMessage().startswith("chembl_request_split")
        for record in caplog.records
    )
