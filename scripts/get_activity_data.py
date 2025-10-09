"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import sys

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from library.cli.entrypoints import activity as _activity  # noqa: E402
from library.pipelines.activity import runner as _activity_runner  # noqa: E402

_MAX_ACTIVITY_ATTRS: dict[str, object] = {
    "MAX_ACTIVITY_CHUNK_SIZE": _activity_runner.MAX_ACTIVITY_CHUNK_SIZE,
    "register_activity_pipeline_hooks": _activity_runner.register_activity_pipeline_hooks,
}

for name, value in _MAX_ACTIVITY_ATTRS.items():
    setattr(_activity, name, value)

_existing_all = tuple(getattr(_activity, "__all__", ()))
_missing_exports = tuple(name for name in _MAX_ACTIVITY_ATTRS if name not in _existing_all)
if _missing_exports:
    __all__ = _existing_all + _missing_exports
    setattr(_activity, "__all__", __all__)
else:
    __all__ = _existing_all

def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.entrypoints.activity`."""

    try:
        return getattr(_activity, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())
