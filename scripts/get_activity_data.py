"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import sys
from importlib import import_module
from types import ModuleType


# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli


def _augment_activity_module(module: ModuleType) -> tuple[str, ...]:
    """Install backwards compatible exports on the activity CLI module."""

    from library.pipelines.activity import runner as activity_runner

    extra_exports: dict[str, object] = {
        "MAX_ACTIVITY_CHUNK_SIZE": activity_runner.MAX_ACTIVITY_CHUNK_SIZE,
        "register_activity_pipeline_hooks": activity_runner.register_activity_pipeline_hooks,
    }

    for name, value in extra_exports.items():
        setattr(module, name, value)

    existing_all = tuple(getattr(module, "__all__", ()))
    missing = tuple(name for name in extra_exports if name not in existing_all)
    if missing:
        module.__all__ = existing_all + missing  # type: ignore[attr-defined]
    return tuple(getattr(module, "__all__", ()))


def _load_activity_module() -> ModuleType:
    module = import_module("library.cli.entrypoints.activity")
    _augment_activity_module(module)
    return module


_activity = _load_activity_module()
__all__ = tuple(getattr(_activity, "__all__", ()))


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
