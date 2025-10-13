"""Unit tests for dictionary path resolution in activity post-processing."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from config.paths import DICTIONARY_DIR
from library.postprocessing import activity_extended
from library.postprocessing.common.runtime import (
    configure_runtime_paths_from_config,
    get_default_export_root,
    override_default_export_root,
)


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


@pytest.mark.unit
def test_current_default_search_dir__prefers_runtime_override(tmp_path: Path) -> None:
    """Runtime export root should override the legacy fallback directory."""

    override_dir = tmp_path / "exports"
    override_dir.mkdir()

    with override_default_export_root(override_dir):
        resolved = activity_extended._current_default_search_dir()
        assert resolved == override_dir

    # Ensure the runtime override has been restored to its previous value.
    assert get_default_export_root() is None


@pytest.mark.unit
def test_configure_runtime_paths_from_config__uses_config_output(
    tmp_path: Path,
) -> None:
    """Configuration objects should populate the runtime export root."""

    export_dir = tmp_path / "configured"
    cfg = SimpleNamespace(io=SimpleNamespace(output_dir=export_dir))

    with override_default_export_root(None):
        configure_runtime_paths_from_config(cfg)
        resolved = activity_extended._current_default_search_dir()
        assert resolved == export_dir
