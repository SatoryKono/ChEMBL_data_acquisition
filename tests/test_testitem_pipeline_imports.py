"""Tests covering import helpers for :mod:`library.testitem_pipeline`."""

from __future__ import annotations

import importlib
from pathlib import Path
import sys

import pytest


def test_load_local_module_missing_file() -> None:
    """Provide a clear error when the optional catalog module is absent."""

    module = importlib.import_module("library.testitem_pipeline")
    missing_name = "missing_catalog"
    qualified = f"{module.__name__}.{missing_name}"
    module_base = Path(module.__file__).with_name(missing_name)
    assert not (module_base.with_suffix(".py").exists())
    assert not (module_base.with_suffix(".pyc").exists())

    sys_modules_before = set(sys.modules)

    with pytest.raises(ModuleNotFoundError) as excinfo:
        module._load_local_module(missing_name)

    message = str(excinfo.value)
    assert qualified in message
    assert f"{missing_name}.py" in message
    assert f"{missing_name}.pyc" in message
    assert qualified not in sys.modules
    assert sys_modules_before <= set(sys.modules)
