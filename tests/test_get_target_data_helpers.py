"""Tests for helper utilities in :mod:`get_target_data`."""

from scripts.get_target_data import _first_token, _pipe_merge


def test_pipe_merge_deduplicates_and_sorts() -> None:
    values = ["b|a", "C|a", None, "", "b"]
    assert _pipe_merge(values) == "C|a|b"


def test_pipe_merge_handles_empty() -> None:
    assert _pipe_merge([None, "", " "]) == ""


def test_first_token_extracts_head() -> None:
    assert _first_token("a|b|c") == "a"
    assert _first_token(None) == ""
    assert _first_token("") == ""
