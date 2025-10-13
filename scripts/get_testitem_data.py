"""Compatibility wrapper exposing :mod:`library.cli.commands.get_testitem_data`."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from importlib import import_module
from pathlib import Path
from types import ModuleType

try:  # pragma: no cover - exercised via CLI invocation
    from scripts._bootstrap import ensure_project_root
except ModuleNotFoundError as exc:  # pragma: no cover - import guard
    if exc.name not in {"scripts", "scripts._bootstrap"}:
        raise
    current_dir = Path(__file__).resolve().parent
    project_root = current_dir.parent
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)
    sys.modules.pop("scripts", None)
    from scripts._bootstrap import ensure_project_root


def _export_module_api(
    module: ModuleType, *, extra: Iterable[str] = ()
) -> tuple[str, ...]:
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


def _load_module() -> ModuleType:
    ensure_project_root(__file__)

    try:
        return import_module("library.cli.commands.get_testitem_data")
    except ModuleNotFoundError as exc:
        if exc.name != "library":
            raise
        ensure_project_root(__file__)
        return import_module("library.cli.commands.get_testitem_data")


_MODULE = _load_module()
__all__ = _export_module_api(_MODULE)


def __getattr__(name: str) -> object:  # pragma: no cover - passthrough helper
    try:
        return getattr(_MODULE, name)
    except AttributeError as exc:  # pragma: no cover - propagate missing attrs
        raise AttributeError(name) from exc


def __dir__() -> list[str]:  # pragma: no cover - passthrough helper
    return sorted({*globals().keys(), *dir(_MODULE)})


class _Adapter(ModuleType):
    """Proxy module syncing attribute writes back to ``_MODULE``."""

    def __getattr__(self, name: str) -> object:  # pragma: no cover - passthrough helper
        return __getattr__(name)

    def __setattr__(
        self, name: str, value: object
    ) -> None:  # pragma: no cover - passthrough helper
        setattr(_MODULE, name, value)
        super().__setattr__(name, value)

    def __dir__(self) -> list[str]:  # pragma: no cover - passthrough helper
        return __dir__()


sys.modules[__name__].__class__ = _Adapter


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_MODULE.main())
