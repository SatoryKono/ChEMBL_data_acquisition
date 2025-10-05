from __future__ import annotations

from collections import deque
from copy import deepcopy

import pandas as pd
import pytest

from library.config import ApiCfg
from library.pipelines.tissue.chembl import TISSUE_BASE_COLUMNS, get_tissues


class _FailingClient:
    """Client stub raising if :meth:`request_json` is invoked."""

    def __init__(self) -> None:
        self.calls: list[str] = []

    def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None = None) -> dict:
        self.calls.append(url)
        raise AssertionError("request_json should not be called")


class _SequenceClient:
    """Client stub returning predetermined payloads in call order."""

    def __init__(self, responses: list[dict[str, object]]) -> None:
        self._responses: deque[dict[str, object]] = deque(deepcopy(responses))
        self.calls: list[dict[str, object]] = []

    def request_json(self, url: str, *, cfg: ApiCfg, timeout: float | None = None) -> dict:
        self.calls.append({"url": url, "timeout": timeout})
        if not self._responses:
            raise AssertionError(f"Unexpected request for {url}")
        return deepcopy(self._responses.popleft())


@pytest.mark.unit
def test_get_tissues__skips_invalid_identifiers() -> None:
    cfg = ApiCfg()
    client = _FailingClient()

    result = get_tissues(["", "#N/A"], cfg=cfg, client=client)

    assert result.empty
    assert list(result.columns) == TISSUE_BASE_COLUMNS
    assert client.calls == []


@pytest.mark.unit
def test_get_tissues__chunks_requests_and_normalises_records() -> None:
    cfg = ApiCfg()
    base = cfg.chembl_base.rstrip("/")
    responses: list[dict[str, object]] = [
        {
            "tissues": [
                {
                    "tissue_chembl_id": "CHEMBL613507",
                    "pref_name": "Liver",
                    "uberon_id": "UBERON:0002107",
                    "efo_id": "EFO:0002067",
                    "bto_id": None,
                    "caloha_id": "TS-0565",
                    "extra": "ignored",
                }
            ],
            "page_meta": {
                "next": (
                    f"{base}/tissue.json?format=json"
                    "&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2&offset=1"
                )
            },
        },
        {
            "tissues": [
                {
                    "tissue_chembl_id": "CHEMBL2109249",
                    "pref_name": "Cerebral cortex",
                    "uberon_id": "UBERON:0000956",
                    "bto_id": "BTO:0000142",
                }
            ],
            "page_meta": {},
        },
        {
            "tissues": [],
            "page_meta": {},
        },
    ]
    client = _SequenceClient(responses)

    frame = get_tissues(
        ["CHEMBL613507", "CHEMBL2109249", "CHEMBL999999"],
        cfg=cfg,
        client=client,
        chunk_size=2,
        timeout=12.5,
    )

    assert list(frame.columns) == TISSUE_BASE_COLUMNS
    assert list(frame["tissue_chembl_id"]) == ["CHEMBL613507", "CHEMBL2109249"]
    assert frame.loc[0, "pref_name"] == "Liver"
    assert frame.loc[1, "pref_name"] == "Cerebral cortex"
    assert pd.isna(frame.loc[0, "bto_id"])
    assert pd.isna(frame.loc[1, "efo_id"])

    urls = [call["url"] for call in client.calls]
    assert urls == [
        (
            f"{base}/tissue.json?format=json"
            "&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2"
        ),
        (
            f"{base}/tissue.json?format=json"
            "&tissue_chembl_id__in=CHEMBL613507,CHEMBL2109249&limit=2&offset=1"
        ),
        (f"{base}/tissue.json?format=json&tissue_chembl_id__in=CHEMBL999999&limit=1"),
    ]
    timeouts = [call["timeout"] for call in client.calls]
    assert timeouts == [12.5, 12.5, 12.5]
