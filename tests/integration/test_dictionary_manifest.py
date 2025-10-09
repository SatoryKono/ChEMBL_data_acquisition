"""Integration smoke tests for the bundled dictionary manifest."""

from __future__ import annotations

import pytest

from library.resources import dictionaries


@pytest.mark.integration
def test_dictionary_manifest__bundled_resources_validate() -> None:
    """Ensure the checked-in manifest matches the bundled resources."""

    dictionaries._env_checksum_allowlist.cache_clear()
    try:
        resources = dictionaries.list_resources()
    finally:
        dictionaries._env_checksum_allowlist.cache_clear()

    assert "dictionary_root" in resources
    for resource in resources.values():
        assert resource.path.exists(), f"missing resource payload: {resource.name}"
        computed = dictionaries._compute_sha256(resource.path)
        assert (
            computed == resource.sha256
        ), f"checksum mismatch for {resource.name}: {computed}"
