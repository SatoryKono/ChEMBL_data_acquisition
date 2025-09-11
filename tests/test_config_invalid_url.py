"""Test configuration validation for malformed URLs."""

from pathlib import Path

import pytest

from library.config import load_config


def test_env_alias_invalid_url(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Invalid base URL supplied via env alias should raise ``ValueError``."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_BASE", "https://")
    with pytest.raises(ValueError, match="api.chembl_base"):
        load_config(path)
