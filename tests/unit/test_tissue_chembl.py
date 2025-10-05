"""Unit tests for :mod:`library.pipelines.tissue.chembl`."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from library.pipelines.tissue import chembl


@dataclass
class _CapturedRequest:
    url: str
    timeout: float | None


class _FakeClient:
    """Minimal stub emulating :class:`ChemblClient`."""

    def __init__(self) -> None:
        self.requests: list[_CapturedRequest] = []

    def request_json(self, url: str, *, cfg, timeout: float | None) -> dict[str, object]:
        self.requests.append(_CapturedRequest(url, timeout))
        if "CHEMBLT1" in url:
            return {
                "tissues": [
                    {
                        "tissue_chembl_id": "CHEMBLT1",
                        "pref_name": "Alpha",
                        "uberon_id": "UBERON:0001",
                        "extra": "ignored",
                    }
                ],
                "page_meta": {"next": None},
            }
        if "CHEMBLT2" in url:
            return {
                "tissues": [
                    {
                        "tissue_chembl_id": "CHEMBLT2",
                        "pref_name": "Beta",
                        "efo_id": "EFO:0001",
                    }
                ],
            }
        return {"tissues": [], "page_meta": {"next": None}}


def test_get_tissues__chunks_and_normalises(cfg) -> None:
    """Multiple identifiers are fetched in chunks and normalised."""

    client = _FakeClient()

    frame = chembl.get_tissues(
        ["CHEMBLT1", "CHEMBLT2"],
        cfg=cfg.api,
        client=client,
        chunk_size=1,
        timeout=12.5,
    )

    assert [req.url for req in client.requests] == [
        f"{cfg.api.chembl_base.rstrip('/')}/tissue.json?format=json&"
        "tissue_chembl_id__in=CHEMBLT1&limit=1",
        f"{cfg.api.chembl_base.rstrip('/')}/tissue.json?format=json&"
        "tissue_chembl_id__in=CHEMBLT2&limit=1",
    ]
    assert all(req.timeout == 12.5 for req in client.requests)
    assert list(frame.columns) == chembl.TISSUE_BASE_COLUMNS
    expected = pd.DataFrame(
        {
            "tissue_chembl_id": ["CHEMBLT1", "CHEMBLT2"],
            "pref_name": ["Alpha", "Beta"],
            "uberon_id": ["UBERON:0001", pd.NA],
            "efo_id": [pd.NA, "EFO:0001"],
            "bto_id": [pd.NA, pd.NA],
            "caloha_id": [pd.NA, pd.NA],
        }
    ).astype(pd.StringDtype())
    assert frame.astype(pd.StringDtype()).equals(expected)

