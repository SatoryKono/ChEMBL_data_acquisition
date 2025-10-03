"""Tests for configuration file discovery helpers."""

from __future__ import annotations

from library.utils import config


def test_resolve_config_path_prefers_repository_default(monkeypatch, tmp_path):
    """Bare config names fall back to the repository config directory."""

    # Simulate running outside the repository tree where ``config.yaml`` is
    # not present in the working directory.
    monkeypatch.chdir(tmp_path)

    resolved = config.resolve_config_path("config.yaml")

    assert resolved == config.DEFAULT_CONFIG_PATH.resolve()
    assert resolved.exists(), "Default configuration must be available"
