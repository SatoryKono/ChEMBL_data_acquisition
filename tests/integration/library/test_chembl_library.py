"""Tests for :mod:`library.integration.chembl_library`."""

from __future__ import annotations

import pytest

from library.integration import chembl_library as cl


def test_chunked_splits_list() -> None:
    items = ["a", "b", "c", "d", "e"]
    chunks = list(cl._chunked(items, 2))
    assert chunks == [["a", "b"], ["c", "d"], ["e"]]


def test_chunked_size_validation() -> None:
    with pytest.raises(ValueError):
        list(cl._chunked(["a"], 0))


def test_chunked_supports_generators() -> None:
    items = (str(i) for i in range(5))
    chunks = list(cl._chunked(items, 2))
    assert chunks == [["0", "1"], ["2", "3"], ["4"]]
