"""Unit tests for :mod:`library.postprocessing.assay_extended`."""

from __future__ import annotations

from pathlib import Path

import pytest

from config.paths import DICTIONARY_DIR
from library.postprocessing import assay_extended


@pytest.mark.unit
def test_resolve_dictionary_root__defaults_to_config_directory(monkeypatch, tmp_path: Path) -> None:
    """Default dictionary path should always be anchored to the config package."""

    # Changing the working directory mimics invocation from an arbitrary location.
    monkeypatch.chdir(tmp_path)

    resolved = assay_extended._resolve_dictionary_root(None)

    assert resolved == DICTIONARY_DIR
    assert resolved.is_absolute()


@pytest.mark.unit
def test_resolve_dictionary_root__accepts_string_paths(tmp_path: Path) -> None:
    """Explicit dictionary_dir arguments should be normalised to ``Path`` instances."""

    provided = tmp_path / "dictionary"
    provided.mkdir()

    resolved = assay_extended._resolve_dictionary_root(str(provided))

    assert resolved == provided
    assert isinstance(resolved, Path)
