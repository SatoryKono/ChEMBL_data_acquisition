"""Unit tests for dictionary path resolution in activity post-processing."""

from __future__ import annotations

from pathlib import Path

import pytest

from config.paths import DICTIONARY_DIR
from library.postprocessing import activity_extended


@pytest.mark.unit
def test_resolve_dictionary_root__defaults_to_config_directory(
    monkeypatch, tmp_path: Path
) -> None:
    """Default dictionary root should always match the config package resource."""

    monkeypatch.chdir(tmp_path)

    resolved = activity_extended._resolve_dictionary_root(None)

    assert resolved == DICTIONARY_DIR
    assert resolved.is_absolute()


@pytest.mark.unit
def test_resolve_dictionary_root__accepts_string_paths(tmp_path: Path) -> None:
    """``dictionary_dir`` parameters passed as strings must resolve correctly."""

    provided = tmp_path / "dictionary"
    provided.mkdir()

    resolved = activity_extended._resolve_dictionary_root(str(provided))

    assert resolved == provided
    assert isinstance(resolved, Path)
