import pytest

from chembl_da.library.config import load_config


def test_load_defaults(tmp_path):
    cfg = load_config(tmp_path / "missing.yaml")
    assert cfg.api.rps == 5
    assert cfg.jobs.concurrency == 8


def test_load_yaml(tmp_path):
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text("api:\n  rps: 2\n")
    cfg = load_config(cfg_file)
    assert cfg.api.rps == 2


def test_env_override(monkeypatch, tmp_path):
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text("api:\n  rps: 2\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "3")
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", "xyz")
    cfg = load_config(cfg_file)
    assert cfg.api.rps == 3
    assert cfg.io.output_dir == "xyz"


def test_validation_error(tmp_path):
    cfg_file = tmp_path / "bad.yaml"
    cfg_file.write_text("api:\n  rps: 0\n")
    with pytest.raises(ValueError):
        load_config(cfg_file)
