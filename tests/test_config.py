from pathlib import Path

import pytest

from library.config import load_config


def test_load_minimal_config(tmp_path: Path) -> None:
    cfg = load_config(tmp_path / "missing.yaml")
    assert cfg.api.rps == 5

    assert cfg.openalex.rps == 4


def test_env_overrides_yaml(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 1\nopenalex:\n  rps: 1\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "3")
    monkeypatch.setenv("CHEMBL_DA__OPENALEX__RPS", "6")
    cfg = load_config(path)
    assert cfg.api.rps == 3
    assert cfg.openalex.rps == 6


def test_alias_env_overrides(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure shorthand environment variables map to the correct fields."""

    path = tmp_path / "cfg.yaml"
    path.write_text(
        "api:\n"
        "  chembl_base: https://www.ebi.ac.uk/chembl/api/data\n"
        "  timeout_connect: 1\n"
        "  timeout_read: 1\n"
    )

    monkeypatch.setenv("CHEMBL_DA_BASE", "https://example.org")
    monkeypatch.setenv("CHEMBL_DA_TIMEOUT_CONNECT", "10")
    monkeypatch.setenv("CHEMBL_DA_TIMEOUT_READ", "20")
    cfg = load_config(path)

    assert cfg.api.chembl_base == "https://example.org"
    assert cfg.api.timeout_connect == 10
    assert cfg.api.timeout_read == 20
    assert cfg.api.rps == 5


def test_cli_overrides_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 1\n")

    monkeypatch.setenv("CHEMBL_DA__API__RPS", "2")
    cfg = load_config(path, cli_overrides={"api.rps": 4})
    assert cfg.api.rps == 4


def test_type_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("api:\n  rps: fast\n")
    with pytest.raises(TypeError, match="api.rps must be int"):
        load_config(path)
