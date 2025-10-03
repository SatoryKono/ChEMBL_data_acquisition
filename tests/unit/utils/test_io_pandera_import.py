from __future__ import annotations

import builtins
import importlib
import sys
from pathlib import Path

import pytest

from library.config import IoCfg


def test_import_without_pandera(monkeypatch: pytest.MonkeyPatch) -> None:
    """Import :mod:`library.io` when ``pandera`` is missing and validate error on schema use."""
    module_name = "library.io"
    orig_import = builtins.__import__

    def fake_import(name: str, *args: object, **kwargs: object):
        if name.startswith("pandera"):
            raise ImportError("boom")
        return orig_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    sys.modules.pop(module_name, None)
    sys.modules.pop("pandera", None)

    io_mod = importlib.import_module(module_name)
    assert io_mod.pa is None

    cfg = IoCfg()
    path = Path("tests/data/pmids.csv")
    with pytest.raises(RuntimeError, match="pandera is required"):
        io_mod.read_csv(path, cfg=cfg, schema=object())

    # restore real module for other tests
    monkeypatch.setattr(builtins, "__import__", orig_import, raising=False)
    sys.modules.pop(module_name, None)
    importlib.import_module(module_name)
