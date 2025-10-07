"""Unit tests for :mod:`library.resources.dictionaries`."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

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
def test_compute_sha256__ignores_os_metadata_case_insensitively(tmp_path: Path) -> None:
    """Common OS metadata files must not impact directory checksums."""

    root = tmp_path / "dictionary"
    root.mkdir()
    (root / "data.csv").write_text("id\n1\n", encoding="utf-8")

    expected = dictionaries._compute_sha256(root)

    metadata_files = [
        "Thumbs.db",
        "thumbs.db",
        "DESKTOP.INI",
        "EhThumbs.DB",
    ]
    for name in metadata_files:
        (root / name).write_bytes(b"\x00\x01")

    (root / "MODULE.PYC").write_bytes(b"\x02\x03")
    pycache = root / "SubDir" / "__PyCache__"
    pycache.mkdir(parents=True)
    (pycache / "ignored.pyc").write_bytes(b"\x04\x05")

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


@pytest.mark.unit
def test_manifest_allows_windows_textmode_checksum() -> None:
    """Ensure the dictionary manifest accepts all known Windows hash variants."""

    manifest_path = Path("config/dictionary/manifest.yaml")
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))

    sha256_values = set(manifest["resources"]["dictionary_root"]["sha256"])

    expected = {
        "efc69f6bb252d68bc7fde11ba98b09b24b0b8fd868fcd6d945eaca76b636f43a",
        "ac67acf2dcd801ffbe9d6e3aa95189af7c3e991fb3ddaaf8aab0be988d7d3224",
        dictionaries.WINDOWS_SPARSE_INDEX_CHECKSUM,
    }

    assert expected.issubset(sha256_values)


@pytest.mark.unit
def test_known_checksum_variants__includes_sparse_index_checksum(tmp_path: Path) -> None:
    """The runtime allow-list should accept the new sparse-index checksum."""

    variants = dictionaries._iter_additional_checksums("dictionary_root", base_dir=tmp_path)

    assert dictionaries.WINDOWS_SPARSE_INDEX_CHECKSUM in variants
