"""Tests for CLI helpers."""
from __future__ import annotations

from pathlib import Path

from main import load_records


def test_load_records_encoding(tmp_path: Path) -> None:
    path = tmp_path / "cp.json"
    content = '{"ids": {}, "pts": {}, "mesh_desc": {}, "mesh_qual": {}}\n'
    path.write_bytes(content.encode("cp1252"))
    data = load_records(path, "cp1252")
    assert data[0]["ids"] == {}
