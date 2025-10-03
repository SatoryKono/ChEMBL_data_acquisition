"""Tests for missing input diagnostics in :mod:`scripts.get_document_data`."""

from __future__ import annotations

from pathlib import Path

from scripts import get_document_data as gdd


def test_missing_input_context_includes_cli_hint(tmp_path: Path) -> None:
    """Ensure the missing input context provides actionable CLI guidance."""

    parent = tmp_path / "inputs"
    parent.mkdir()
    candidate = parent / "document1.csv"
    candidate.write_text("", encoding="utf-8")

    context = gdd._build_missing_input_context(parent / "document.csv")

    assert context["did_you_mean"] == str(candidate)
    assert context["cli_hint"] == f"--input \"{candidate}\""
    assert context["suggestions"] == [str(candidate)]
