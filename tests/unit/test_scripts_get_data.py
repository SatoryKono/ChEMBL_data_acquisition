from __future__ import annotations

import importlib
import runpy
import sys

import pytest


@pytest.mark.unit
def test_wrapper_main__aliases_library_main():
    """Wrapper should expose the library CLI implementation transparently."""

    sys.modules.pop("scripts.get_data", None)

    wrapper = importlib.import_module("scripts.get_data")
    library_module = importlib.import_module("library.cli.commands.get_data")

    assert wrapper.main is library_module.main


@pytest.mark.unit
def test_main_entrypoint__delegates_to_library_main(monkeypatch):
    """Executing the wrapper as a script should invoke the library main."""

    library_module = importlib.import_module("library.cli.commands.get_data")

    calls: dict[str, object] = {}

    def _fake_main(argv=None):
        calls["argv"] = argv
        return 17

    monkeypatch.setattr(library_module, "main", _fake_main)
    sys.modules.pop("scripts.get_data", None)

    previous_main = sys.modules.get("__main__")
    try:
        with pytest.raises(SystemExit) as excinfo:
            runpy.run_module("scripts.get_data", run_name="__main__")
    finally:
        if previous_main is not None:
            sys.modules["__main__"] = previous_main
        else:
            sys.modules.pop("__main__", None)

    assert excinfo.value.code == 17
    assert calls["argv"] is None
