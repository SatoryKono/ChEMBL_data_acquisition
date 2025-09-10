from __future__ import annotations

from pathlib import Path

import pytest

from library.config import load_config


def test_defaults_without_file(tmp_path: Path) -> None:
    cfg = load_config(tmp_path / "missing.yaml")
    assert cfg.api.rps == 5
    assert cfg.io.output_dir == "data/output"


def test_yaml_overrides(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 2\n")
    cfg = load_config(path)
    assert cfg.api.rps == 2


def test_env_overrides(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 1\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "7")
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", "out")
    cfg = load_config(path)
    assert cfg.api.rps == 7
    assert cfg.io.output_dir == "out"


def test_init_paths(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("init:\n  same_doc: custom_same.xlsx\n  all_doc: custom_all.xlsx\n")
    monkeypatch.setenv("CHEMBL_DA__INIT__ALL_DOC", "env_all.xlsx")
    cfg = load_config(path)
    assert cfg.init.same_doc == "custom_same.xlsx"
    assert cfg.init.all_doc == "env_all.xlsx"


def test_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("api:\n  rps: 0\n")
    with pytest.raises(ValueError):
        load_config(path)
