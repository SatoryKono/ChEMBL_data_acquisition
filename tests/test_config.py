from pathlib import Path

import pytest

from library.config import ensure_dirs, load_config


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


def test_cli_path_override(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    out = tmp_path / "out"
    cfg = load_config(path, cli_overrides={"io.output_dir": str(out)})
    assert cfg.io.output_dir == out


def test_type_validation(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("api:\n  rps: fast\n")
    with pytest.raises(TypeError, match="api.rps must be int"):
        load_config(path)


def test_missing_dirs_raise(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(tmp_path / "out"))
    monkeypatch.setenv("CHEMBL_DA__IO__CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("CHEMBL_DA__IO__EXIST_OK", "false")
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "cfg.yaml")


def test_ensure_dirs_creates(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "out"
    cache = tmp_path / "cache"
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(out))
    monkeypatch.setenv("CHEMBL_DA__IO__CACHE_DIR", str(cache))
    cfg = load_config(tmp_path / "cfg.yaml")
    assert not out.exists() and not cache.exists()
    ensure_dirs(cfg)
    assert out.is_dir() and cache.is_dir()
