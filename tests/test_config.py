"""Tests for the configuration utilities."""

from pathlib import Path

import pytest

from library.config import (
    APISettings,
    Config,
    ConfigError,
    OutputPaths,
    TimeoutSettings,
    load_config,
)


DATA_DIR = Path(__file__).resolve().parent / "data"


def test_load_config_reads_yaml() -> None:
    path = DATA_DIR / "sample_config.yaml"
    cfg = load_config(path)
    assert cfg.api.chembl_base_url == "https://example.org/chembl"
    assert cfg.timeouts.connect == 5


def test_load_config_missing_returns_default(tmp_path: Path) -> None:
    path = tmp_path / "missing.yaml"
    cfg = load_config(path)
    assert cfg == Config()


def test_negative_timeout_raises() -> None:
    with pytest.raises(ConfigError):
        Config(timeouts=TimeoutSettings(connect=-1, read=10))


def test_invalid_url_raises() -> None:
    api = APISettings(chembl_base_url="not-a-url")
    with pytest.raises(ConfigError):
        Config(api=api)


def test_non_writable_directory_raises(tmp_path: Path) -> None:
    bad_path = tmp_path / "file"
    bad_path.write_text("x")
    with pytest.raises(ConfigError):
        Config(output=OutputPaths(data_dir=bad_path))
