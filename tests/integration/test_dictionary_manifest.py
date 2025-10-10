"""Integration smoke tests for the bundled dictionary manifest."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from config.paths import DICTIONARY_DIR
from library.resources import dictionaries


def _copy_dictionary_bundle(tmp_path: Path) -> Path:
    destination = tmp_path / "dictionary"
    shutil.copytree(DICTIONARY_DIR, destination)
    return destination


@pytest.mark.integration
@pytest.mark.smoke
def test_dictionary_manifest__bundled_resources_validate(tmp_path: Path) -> None:
    """Ensure the checked-in manifest matches the bundled resources."""

    dictionary_root = _copy_dictionary_bundle(tmp_path)
    dictionaries._env_checksum_allowlist.cache_clear()
    try:
        resources = dictionaries.list_resources(base_dir=dictionary_root)
    finally:
        dictionaries._env_checksum_allowlist.cache_clear()

    assert "dictionary_root" in resources
    for resource in resources.values():
        assert resource.path.exists(), f"missing resource payload: {resource.name}"
        computed = dictionaries._compute_sha256(resource.path)
        assert (
            computed == resource.sha256
        ), f"checksum mismatch for {resource.name}: {computed}"
