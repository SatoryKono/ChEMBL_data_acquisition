"""Unit tests for :func:`library.pipelines.tissue.chembl.get_tissues`."""

from __future__ import annotations

from collections import deque
from typing import Any
from unittest.mock import create_autospec
from urllib.parse import parse_qs, urlsplit

import pandas as pd
import pytest

from library.clients import ChemblClient
from library.config import ApiCfg
from library.pipelines.tissue.chembl import TISSUE_BASE_COLUMNS, get_tissues


def _make_payload(
    *,
    tissues: list[dict[str, Any]],
    next_token: str | None = None,
) -> dict[str, Any]:
    """Return a mock payload with optional pagination metadata."""

    page_meta: dict[str, Any] = {}
    if next_token is not None:
        page_meta["next"] = next_token
    return {"tissues": tissues, "page_meta": page_meta}


@pytest.mark.unit
def test_get_tissues__paginates_and_preserves_order() -> None:
    """The client is called with chunked URLs and honours custom timeouts."""

    cfg = ApiCfg(chembl_base="https://example.test/custom", timeout_read=12.5)
    client = create_autospec(ChemblClient, instance=True)
    client.request_json.side_effect = deque(
        [
            _make_payload(
                tissues=[
                    {
                        "tissue_chembl_id": "TIS_1",
                        "pref_name": "Alpha",
                        "uberon_id": "UBERON:1",
                    }
                ],
                next_token="/tissue.json?format=json&offset=1",
            ),
            _make_payload(
                tissues=[
                    {
                        "tissue_chembl_id": "TIS_2",
                        "pref_name": "Beta",
                        "uberon_id": "UBERON:2",
                    }
                ]
            ),
            _make_payload(
                tissues=[
                    {
                        "tissue_chembl_id": "TIS_3",
                        "pref_name": "Gamma",
                        "uberon_id": "UBERON:3",
                    },
                    {
                        "tissue_chembl_id": "TIS_4",
                        "pref_name": "Delta",
                        "uberon_id": "UBERON:4",
                    },
                ]
            ),
        ]
    )

    df = get_tissues(
        ["TIS_1", "TIS_2", "TIS_3", "TIS_4"],
        cfg=cfg,
        client=client,
        chunk_size=2,
        timeout=7.0,
    )

    urls = [call.args[0] for call in client.request_json.call_args_list]
    timeouts = [call.kwargs["timeout"] for call in client.request_json.call_args_list]
    cfgs = [call.kwargs["cfg"] for call in client.request_json.call_args_list]

    assert urls == [
        "https://example.test/custom/tissue.json?format=json"
        "&tissue_chembl_id__in=TIS_1,TIS_2&limit=2",
        "https://example.test/tissue.json?format=json&offset=1",
        "https://example.test/custom/tissue.json?format=json"
        "&tissue_chembl_id__in=TIS_3,TIS_4&limit=2",
    ]
    assert timeouts == [7.0, 7.0, 7.0]
    assert all(obj is cfg for obj in cfgs)
    assert df["tissue_chembl_id"].tolist() == ["TIS_1", "TIS_2", "TIS_3", "TIS_4"]


@pytest.mark.unit
def test_get_tissues__normalises_missing_columns() -> None:
    """Missing keys from the payload are filled with ``pd.NA``."""

    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = create_autospec(ChemblClient, instance=True)
    client.request_json.side_effect = [
        _make_payload(
            tissues=[
                {
                    "tissue_chembl_id": "TIS_MISSING",
                    "pref_name": None,
                    "uberon_id": "",
                    "caloha_id": None,
                }
            ]
        )
    ]

    df = get_tissues(
        ["TIS_MISSING"],
        cfg=cfg,
        client=client,
        chunk_size=5,
    )

    assert df.columns.tolist() == TISSUE_BASE_COLUMNS
    row = df.iloc[0]
    assert row["tissue_chembl_id"] == "TIS_MISSING"
    assert pd.isna(row["pref_name"])
    assert row["uberon_id"] == ""
    assert pd.isna(row["efo_id"])
    assert pd.isna(row["bto_id"])
    assert pd.isna(row["caloha_id"])


@pytest.mark.unit
def test_get_tissues__handles_empty_or_null_payload() -> None:
    """Empty ``tissues`` payloads yield an empty dataframe with base columns."""

    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = create_autospec(ChemblClient, instance=True)
    client.request_json.side_effect = [
        {"tissues": None, "page_meta": {}},
    ]

    df = get_tissues(
        ["CHEMBLT1"],
        cfg=cfg,
        client=client,
        chunk_size=3,
    )

    client.request_json.assert_called_once()
    assert df.empty
    assert df.columns.tolist() == TISSUE_BASE_COLUMNS


@pytest.mark.unit
def test_get_tissues__skips_invalid_identifiers() -> None:
    """Blank identifiers are ignored without contacting the API."""

    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = create_autospec(ChemblClient, instance=True)

    df = get_tissues(
        ["", "#N/A"],
        cfg=cfg,
        client=client,
        chunk_size=3,
    )

    client.request_json.assert_not_called()
    assert df.empty
    assert df.columns.tolist() == TISSUE_BASE_COLUMNS


@pytest.mark.unit
def test_get_tissues__caps_chunk_size_to_service_limit() -> None:
    """Chunk size above API limits is clamped to the maximum of 1000."""

    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=5.0)
    client = create_autospec(ChemblClient, instance=True)
    client.request_json.return_value = {"tissues": [], "page_meta": {}}

    identifiers = [f"TIS_{index}" for index in range(1050)]
    df = get_tissues(
        identifiers,
        cfg=cfg,
        client=client,
        chunk_size=5000,
        timeout=None,
    )

    assert df.empty
    assert df.columns.tolist() == TISSUE_BASE_COLUMNS

    urls = [call.args[0] for call in client.request_json.call_args_list]
    assert len(urls) == 2
    first_limit = parse_qs(urlsplit(urls[0]).query).get("limit", [])
    second_limit = parse_qs(urlsplit(urls[1]).query).get("limit", [])
    assert first_limit == ["1000"]
    assert second_limit == ["50"]
