from __future__ import annotations

import sys
from pathlib import Path

import pytest

from library.config import load_config


def _write_config(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(text, encoding="utf8")
    return path


def test_load_config_from_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_file = _write_config(tmp_path, "api:\n  timeout: 5\n")
    monkeypatch.delenv("CHEMBL_TIMEOUT", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog"])
    cfg = load_config(config_file)
    assert cfg.api.timeout == 5


def test_env_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_file = _write_config(tmp_path, "api:\n  timeout: 5\n")
    monkeypatch.setenv("CHEMBL_TIMEOUT", "10")
    monkeypatch.setattr(sys, "argv", ["prog"])
    cfg = load_config(config_file)
    assert cfg.api.timeout == 10


def test_cli_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config_file = _write_config(tmp_path, "api:\n  timeout: 5\n")
    monkeypatch.delenv("CHEMBL_TIMEOUT", raising=False)
    monkeypatch.setattr(sys, "argv", ["prog", "--timeout", "15"])
    cfg = load_config(config_file)
    assert cfg.api.timeout == 15
