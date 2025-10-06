"""Tests for the legacy ``library.utils.config`` compatibility layer."""

from __future__ import annotations

import importlib

import pytest


@pytest.mark.unit
def test_utils_config__reexports():
    utils_config = importlib.import_module("library.utils.config")

    from library.config import (
        Config,
        ConfigError,
        ConfigLoaderError,
        ConfigMetadata,
        ConfigSource,
        ensure_dirs,
        load_config,
        load_yaml_config,
        print_config,
        resolve_config_path,
    )
    from library.config.loader import DEFAULT_CONFIG_PATH

    assert utils_config.Config is Config
    assert utils_config.ConfigError is ConfigError
    assert utils_config.ConfigLoaderError is ConfigLoaderError
    assert utils_config.ConfigMetadata is ConfigMetadata
    assert utils_config.ConfigSource is ConfigSource
    assert utils_config.ensure_dirs is ensure_dirs
    assert utils_config.load_config is load_config
    assert utils_config.load_yaml_config is load_yaml_config
    assert utils_config.print_config is print_config
    assert utils_config.resolve_config_path is resolve_config_path
    assert utils_config.DEFAULT_CONFIG_PATH == DEFAULT_CONFIG_PATH
