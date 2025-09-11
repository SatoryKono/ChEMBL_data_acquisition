"""Test configuration validation for malformed URLs."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from library.config import load_config


def test_env_alias_invalid_url(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Invalid base URL supplied via env alias should raise ``ValueError``."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_BASE", "https://")
    monkeypatch.setenv(
        "CHEMBL_DA__API__USER_AGENT", "test-agent/1.0 (mailto:test@example.org)"
    )
    with pytest.raises(ValidationError, match="api.chembl_base"):
        load_config(path)
