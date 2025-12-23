from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chembl_da.library.csv_utils import sha256_file


def test_sha256_file(tmp_path: Path) -> None:
    data = b"hello world"
    file_path = tmp_path / "sample.txt"
    file_path.write_bytes(data)
    expected = hashlib.sha256(data).hexdigest()
    assert sha256_file(file_path) == expected


def test_sha256_file_missing(tmp_path: Path) -> None:
    path = tmp_path / "missing.txt"
    with pytest.raises(FileNotFoundError) as excinfo:
        sha256_file(path)
    assert "file not found" in str(excinfo.value)
    assert str(path) in str(excinfo.value)
