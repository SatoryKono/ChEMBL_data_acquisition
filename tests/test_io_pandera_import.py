from __future__ import annotations

import builtins
import importlib
import sys

import pytest


def test_pandera_import_error(monkeypatch: pytest.MonkeyPatch) -> None:
    """Raise a helpful error when pandera import fails."""
    module_name = "library.io"
    orig_import = builtins.__import__

    def fake_import(name: str, *args: object, **kwargs: object):
        if name == "pandera.pandas":
            raise TypeError("boom")
        return orig_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    sys.modules.pop(module_name, None)

    with pytest.raises(RuntimeError, match="pandera is required"):
        importlib.import_module(module_name)

    # restore real module for other tests
    monkeypatch.setattr(builtins, "__import__", orig_import, raising=False)
    sys.modules.pop(module_name, None)
    importlib.import_module(module_name)
