import sys
from pathlib import Path

import pytest

sys.path.append(str(Path(__file__).resolve().parent.parent))

from mylib.io_utils import read_input


def _write_cp1251(path: Path) -> None:
    data = "record_id,doi\n1,тест"
    path.write_text(data, encoding="cp1251")


def test_auto_detect_encoding(tmp_path: Path) -> None:
    file_path = tmp_path / "cp.csv"
    _write_cp1251(file_path)
    df = read_input(file_path, sep=",", encoding=None)
    assert df.loc[0, "doi"] == "тест"


def test_wrong_encoding_raises(tmp_path: Path) -> None:
    file_path = tmp_path / "cp.csv"
    _write_cp1251(file_path)
    with pytest.raises(ValueError):
        read_input(file_path, sep=",", encoding="utf-8")
