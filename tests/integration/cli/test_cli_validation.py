import pytest

from library.cli import build_parser


def test_chunk_size_positive() -> None:
    parser, _ = build_parser("desc", column="id")
    args = parser.parse_args(["--chunk-size", "3"])
    assert args.chunk_size == 3


def test_chunk_size_non_positive() -> None:
    parser, _ = build_parser("desc", column="id")
    with pytest.raises(SystemExit):
        parser.parse_args(["--chunk-size", "0"])
