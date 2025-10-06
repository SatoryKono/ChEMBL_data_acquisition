from __future__ import annotations

import importlib

import pytest


@pytest.mark.integration
def test_config_public_api__reexports():
    package = importlib.import_module("library.config")
    model = importlib.import_module("library.config.models")
    loader = importlib.import_module("library.config.loader")
    runtime = importlib.import_module("library.config.runtime")

    assert package.Config is model.Config
    assert package.ConfigMetadata is model.ConfigMetadata
    assert package.load_config is loader.load_config
    assert package.ensure_dirs is loader.ensure_dirs
    assert package.session_with_retry is runtime.session_with_retry
