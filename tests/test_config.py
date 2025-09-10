from pathlib import Path

import pytest

from library.config import ConfigError, ensure_dirs, load_config


def test_missing_config_raises(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="configuration file not found"):
        load_config(tmp_path / "missing.yaml")


def test_env_hierarchical_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    out = tmp_path / "out"
    monkeypatch.setenv("CHEMBL_DA__IO__OUTPUT_DIR", str(out))
    cfg = load_config(path)
    assert cfg.io.output_dir == out


def test_short_alias_rps(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA_RPS", "9")
    cfg = load_config(path)
    assert cfg.api.rps == 9


def test_cli_override(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path, cli_overrides={"api.rps": 7})
    assert cfg.api.rps == 7


def test_invalid_rps(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: -5\n")
    with pytest.raises(ConfigError):
        load_config(path)


def test_invalid_timeout_type(tmp_path: Path) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  timeout_read: abc\n")
    with pytest.raises(ConfigError, match="api.timeout_read"):
        load_config(path)


@pytest.mark.parametrize(
    "env_val, expected",
    [("true", True), ("false", False)],
)
def test_env_bool_conversion(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, env_val: str, expected: bool
) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("io:\n  exist_ok: false\n")
    monkeypatch.setenv("CHEMBL_DA__IO__EXIST_OK", env_val)
    cfg = load_config(path)
    assert cfg.io.exist_ok is expected


def test_invalid_bool_env_var(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    monkeypatch.setenv("CHEMBL_DA__IO__EXIST_OK", "maybe")
    with pytest.raises(ConfigError, match="Invalid boolean value"):
        load_config(path)


def test_auto_create_dirs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "out"
    cache = tmp_path / "cache"
    init_out = tmp_path / "init_out"
    monkeypatch.setenv("CHEMBL_DA_OUTDIR", str(out))
    monkeypatch.setenv("CHEMBL_DA__IO__CACHE_DIR", str(cache))
    monkeypatch.setenv("CHEMBL_DA__INIT__OUTPUT_DIR", str(init_out))
    path = tmp_path / "cfg.yaml"
    path.write_text("")
    cfg = load_config(path)
    assert out.is_dir() and cache.is_dir() and init_out.is_dir()
    ensure_dirs(cfg)  # idempotent
    assert out.is_dir() and cache.is_dir() and init_out.is_dir()


def test_priority_order(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = tmp_path / "cfg.yaml"
    path.write_text("api:\n  rps: 10\n")
    monkeypatch.setenv("CHEMBL_DA__API__RPS", "20")
    cfg = load_config(path, cli_overrides={"api.rps": 30})
    assert cfg.api.rps == 30
