"""Test configuration validation for malformed URLs."""

from pathlib import Path

import pytest
from library.config import ConfigError, load_config


def test_env_alias_invalid_url(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Invalid base URL supplied via env alias should raise ``ValueError``."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_BASE", "https://")
    monkeypatch.setenv(
        "CHEMBL_DA__SOURCES__CHEMBL__API__USER_AGENT",
        "test-agent/1.0 (mailto:test@example.org)",
    )
    with pytest.raises(ConfigError, match="api.chembl_base"):
        load_config(path)


def test_openalex_placeholder_mailto_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Placeholder OpenAlex email provided via env should raise ``ConfigError``."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_OPENALEX_MAILTO", "contact@example.org")

    with pytest.raises(ConfigError, match="placeholder domain"):
        load_config(path)


def test_crossref_placeholder_mailto_env(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Placeholder CrossRef email provided via env should raise ``ConfigError``."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_CROSSREF_MAILTO", "contact@example.org")

    with pytest.raises(ConfigError, match="placeholder domain"):
        load_config(path)
