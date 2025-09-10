from pathlib import Path
import logging

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


def test_retry_and_log_aliases(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """New environment variable aliases should override retry and log defaults."""

    path = tmp_path / "cfg.yaml"
    path.write_text("")

    monkeypatch.setenv("CHEMBL_DA_RETRY_MAX_ATTEMPTS", "10")
    monkeypatch.setenv("CHEMBL_DA_RETRY_BACKOFF_FACTOR", "2.0")
    monkeypatch.setenv("CHEMBL_DA_LOG_FORMAT", "%(levelname)s")

    cfg = load_config(path)

    assert cfg.retry.max_attempts == 10
    assert cfg.retry.backoff_factor == 2.0
    assert cfg.log.format == "%(levelname)s"


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


def test_unknown_key_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("unknown: 1\napi:\n  rps: 1\n")
    with caplog.at_level(logging.WARNING, logger="library.config"):
        load_config(path)
    assert "Unknown configuration key" in caplog.text


def test_unknown_key_error(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("unknown: 1\n")
    with pytest.raises(ValueError, match="Unknown configuration key"):
        load_config(path, strict=True)


def test_yaml_error_includes_path(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api: [\n")  # malformed YAML
    with pytest.raises(ValueError) as excinfo:
        load_config(path)
    msg = str(excinfo.value)
    assert str(path) in msg
    assert "while parsing" in msg
