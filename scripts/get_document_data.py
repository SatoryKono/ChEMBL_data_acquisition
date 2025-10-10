"""Compatibility shim for :mod:`library.cli.commands.get_document_data`."""

from __future__ import annotations

import sys
from importlib import import_module

# ruff: noqa: E402  # bootstrap updates sys.path for direct script execution
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a package module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

_LEGACY_NAME = __name__
_TARGET_MODULE_NAME = "library.cli.commands.get_document_data"
_CLI_MODULE = import_module(_TARGET_MODULE_NAME)

# Re-export public and private helpers for legacy callers before swapping the module
# object so existing imports keep functioning.
__doc__ = _CLI_MODULE.__doc__
if hasattr(_CLI_MODULE, "__all__"):
    __all__ = _CLI_MODULE.__all__  # type: ignore[assignment]

for _name, _value in vars(_CLI_MODULE).items():
    if _name.startswith("__"):
        continue
    globals()[_name] = _value

# Make ``scripts.get_document_data`` resolve to the CLI implementation module.
sys.modules[_LEGACY_NAME] = _CLI_MODULE

del _name, _value

if __name__ == "__main__":
    sys.exit(_CLI_MODULE.main())
