"""Compatibility wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from importlib import import_module
from types import ModuleType
from typing import Iterable, TYPE_CHECKING

# ruff: noqa: E402  # bootstrap alters import order for script compatibility
if TYPE_CHECKING:
    from . import _bootstrap as _bootstrap_module
elif __package__ in {None, ""}:
    import _bootstrap as _bootstrap_module  # pragma: no cover - CLI fallback
else:  # pragma: no cover - executed when imported as a module
    from . import _bootstrap as _bootstrap_module

bootstrap_cli = _bootstrap_module.bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli
del _bootstrap_module


def _export_module_api(module: ModuleType, *, extra: Iterable[str] = ()) -> tuple[str, ...]:
    """Expose ``module`` attributes in the wrapper namespace."""

    exported: dict[str, object] = {}
    for name, value in module.__dict__.items():
        if name.startswith("__"):
            continue
        exported[name] = value

    for name in extra:
        if hasattr(module, name):
            exported[name] = getattr(module, name)

    globals().update(exported)

    module_all = tuple(getattr(module, "__all__", ()))
    if module_all:
        ordered = list(module_all)
        for name in exported:
            if name not in ordered:
                ordered.append(name)
    else:
        ordered = list(exported)
    return tuple(ordered)


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
    module = import_module("library.cli.commands.get_activity_data")
    return module


_MODULE = _load_activity_module()
__all__ = _export_module_api(_MODULE)


def __getattr__(name: str) -> object:
    """Proxy attribute access to :mod:`library.cli.commands.get_activity_data`."""

    try:
        return getattr(_MODULE, name)
    except AttributeError as exc:  # pragma: no cover - passthrough for missing attrs
        raise AttributeError(name) from exc


def __dir__() -> list[str]:  # pragma: no cover - passthrough helper
    return sorted({*globals().keys(), *dir(_MODULE)})


class _Adapter(ModuleType):
    """Proxy module syncing attribute writes back to ``_MODULE``."""

    def __getattr__(self, name: str) -> object:  # pragma: no cover - passthrough helper
        return __getattr__(name)

    def __setattr__(self, name: str, value: object) -> None:  # pragma: no cover - passthrough helper
        setattr(_MODULE, name, value)
        super().__setattr__(name, value)

    def __dir__(self) -> list[str]:  # pragma: no cover - passthrough helper
        return __dir__()


sys.modules[__name__].__class__ = _Adapter


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_MODULE.main())
