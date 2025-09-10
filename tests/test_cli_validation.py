import pytest

from pathlib import Path

from library.cli import build_parser
from library.config import DEFAULT_PATH


def test_chunk_size_positive() -> None:
    parser = build_parser("desc", column="id")
    args = parser.parse_args(["--chunk-size", "3"])
    assert args.chunk_size == 3


def test_chunk_size_non_positive() -> None:
    parser = build_parser("desc", column="id")
    with pytest.raises(SystemExit):
        parser.parse_args(["--chunk-size", "0"])


def test_config_default_path() -> None:
    parser = build_parser("desc", column="id")
    args = parser.parse_args([])
    assert args.config == Path(DEFAULT_PATH)
