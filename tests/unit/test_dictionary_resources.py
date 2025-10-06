"""Unit tests for :mod:`library.resources.dictionaries`."""

from __future__ import annotations

from pathlib import Path

import pytest

from library.resources import dictionaries


@pytest.mark.unit
def test_compute_sha256__ignores_crlf_in_files(tmp_path: Path) -> None:
    """CRLF newlines should not affect the checksum of individual files."""

    lf_path = tmp_path / "lf.csv"
    lf_path.write_text("id,name\n1,example\n", encoding="utf-8")

    crlf_path = tmp_path / "crlf.csv"
    crlf_path.write_bytes("id,name\r\n1,example\r\n".encode("utf-8"))

    assert dictionaries._compute_sha256(lf_path) == dictionaries._compute_sha256(crlf_path)


@pytest.mark.unit
def test_compute_sha256__ignores_crlf_in_directories(tmp_path: Path) -> None:
    """Directory checksums should be stable across newline conventions."""

    lf_dir = tmp_path / "lf"
    lf_dir.mkdir()
    (lf_dir / "data.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    (lf_dir / "nested").mkdir()
    (lf_dir / "nested" / "more.csv").write_text("x\ny\n", encoding="utf-8")

    crlf_dir = tmp_path / "crlf"
    crlf_dir.mkdir()
    (crlf_dir / "data.csv").write_bytes("a,b\r\n1,2\r\n".encode("utf-8"))
    (crlf_dir / "nested").mkdir()
    (crlf_dir / "nested" / "more.csv").write_bytes("x\r\ny\r\n".encode("utf-8"))

    assert dictionaries._compute_sha256(lf_dir) == dictionaries._compute_sha256(crlf_dir)


@pytest.mark.unit
def test_compute_sha256__ignores_bytecode_artifacts(tmp_path: Path) -> None:
    """Ensure cached Python bytecode files do not affect resource checksums."""

    root = tmp_path / "dictionary"
    root.mkdir()

    data_file = root / "data.csv"
    data_file.write_text("id\n1\n", encoding="utf-8")

    expected = dictionaries._compute_sha256(root)

    pycache = root / "__pycache__"
    pycache.mkdir()
    (pycache / "module.cpython-313.pyc").write_bytes(b"\x00\x01")
    (pycache / "module.pyo").write_bytes(b"\x02\x03")
    (root / "standalone.pyc").write_bytes(b"\x04\x05")

    assert dictionaries._compute_sha256(root) == expected


@pytest.mark.unit
def test_compute_sha256__normalises_crlf_in_non_utf8_files(tmp_path: Path) -> None:
    """CRLF handling must be deterministic for legacy encodings (e.g. Latin-1)."""

    latin1_bytes_lf = b"id;name\n1;\xb1\n"
    latin1_bytes_crlf = b"id;name\r\n1;\xb1\r\n"

    lf_path = tmp_path / "latin1_lf.csv"
    lf_path.write_bytes(latin1_bytes_lf)

    crlf_path = tmp_path / "latin1_crlf.csv"
    crlf_path.write_bytes(latin1_bytes_crlf)

    assert dictionaries._compute_sha256(lf_path) == dictionaries._compute_sha256(crlf_path)
