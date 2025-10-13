"""Unit tests for :func:`library.pipelines.cellline.chembl.get_cell_lines`."""

from __future__ import annotations

from unittest.mock import create_autospec
from urllib.parse import parse_qs, urlsplit

import pytest

from library.clients import ChemblClient
from library.config import ApiCfg
from library.pipelines.cellline.chembl import CELL_LINE_BASE_COLUMNS, get_cell_lines


@pytest.mark.unit
def test_get_cell_lines__caps_chunk_size_to_service_limit() -> None:
    """Chunk sizes larger than the ChEMBL limit are clamped to 1000."""

    cfg = ApiCfg(chembl_base="https://example.test/base", timeout_read=4.0)
    client = create_autospec(ChemblClient, instance=True)
    client.request_json.return_value = {"cell_lines": [], "page_meta": {}}

    identifiers = [f"CELL_{index}" for index in range(1100)]
    df = get_cell_lines(
        identifiers,
        cfg=cfg,
        client=client,
        chunk_size=4096,
    )

    assert df.empty
    assert df.columns.tolist() == CELL_LINE_BASE_COLUMNS

    urls = [call.args[0] for call in client.request_json.call_args_list]
    assert len(urls) == 2
    first_limit = parse_qs(urlsplit(urls[0]).query).get("limit", [])
    second_limit = parse_qs(urlsplit(urls[1]).query).get("limit", [])
    assert first_limit == ["1000"]
    assert second_limit == ["100"]
