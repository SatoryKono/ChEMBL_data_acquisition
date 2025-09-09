from pathlib import Path

import pytest

from chembl_da.library.config import load_config


def test_load_defaults(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cfg = load_config("nonexistent.yaml")
    assert cfg.api.timeout_connect == 5
    assert cfg.jobs.chunk_size == 500


def test_load_from_yaml(tmp_path: Path):
    yaml_path = tmp_path / "conf.yaml"
    yaml_path.write_text(
        """
api:
  rps: 10
io:
  output_dir: "out"
"""
    )
    cfg = load_config(str(yaml_path))
    assert cfg.api.rps == 10
    assert cfg.io.output_dir == "out"


def test_env_override_precedence(tmp_path: Path, monkeypatch):
    yaml_path = tmp_path / "conf.yaml"
    yaml_path.write_text("api:\n  rps: 8\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "9")
    monkeypatch.setenv("CHEMBL_DA_RPS", "11")
    cfg = load_config(str(yaml_path))
    assert cfg.api.rps == 11


def test_env_cast_int(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("CHEMBL_DA__JOBS__CHUNK_SIZE", "42")
    cfg = load_config(str(tmp_path / "missing.yaml"))
    assert cfg.jobs.chunk_size == 42


def test_validation_negative(tmp_path: Path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("jobs:\n  chunk_size: -1\n")
    with pytest.raises(ValueError):
        load_config(str(bad))
