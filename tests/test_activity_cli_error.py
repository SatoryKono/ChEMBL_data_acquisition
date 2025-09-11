import argparse
import logging
from pathlib import Path

import requests

import get_activity_data as gad
from library import chembl_library as cl, io
from library.config import Config


def test_run_chembl_handles_request_error(monkeypatch, caplog) -> None:
    args = argparse.Namespace(
        input_csv=Path("dummy.csv"),
        output_csv=None,
        column="activity_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
        timeout=30.0,
    )
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: ["1"])

    def _raise(*_a, **_k):
        raise requests.RequestException("boom")

    monkeypatch.setattr(cl, "get_activities", _raise)
    with caplog.at_level(logging.ERROR):
        result = gad.run_chembl(Config(), args)
    assert result == 1
    assert "boom" in caplog.text
