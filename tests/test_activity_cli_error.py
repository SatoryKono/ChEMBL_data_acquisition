import argparse
import json
from pathlib import Path
import io

import requests

from library.cli import LoggerConfig, configure_logger
import get_activity_data as gad
from library import chembl_library as cl, io as lib_io
from library.config import Config


def test_run_chembl_handles_request_error(monkeypatch) -> None:
    args = argparse.Namespace(
        input_csv=Path("dummy.csv"),
        output_csv=None,
        column="activity_id",
        sep=",",
        encoding="utf8",
        chunk_size=5,
        timeout=30.0,
        limit=None,
        dry_run=False,
    )
    monkeypatch.setattr(lib_io, "read_ids", lambda *a, **k: iter(["1"]))

    def _raise(*_a, **_k):
        raise requests.RequestException("boom")

    monkeypatch.setattr(cl, "get_activities", _raise)
    buf = io.StringIO()
    configure_logger(LoggerConfig(stream=buf))
    result = gad.run_chembl(Config(), args)
    assert result == 1
    lines = buf.getvalue().splitlines()
    assert lines
    record = json.loads(lines[-1])
    assert "boom" in record.get("msg", "")
