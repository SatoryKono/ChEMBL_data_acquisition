"""Top-level re-exports for the :mod:`library` package."""

from __future__ import annotations

from importlib import import_module
from typing import Any

_LAZY_SUBMODULES = {
    "cli",
    "cli_utils",
    "io",
    "offline",
    "postprocess",
    "qa",
    "schemas",
    "testitem_pipeline",
    "validation",
}

_EXPORTS: dict[str, str] = {
    "SidecarErrors": "library.sidecar:SidecarErrors",
    "resolve_failure_chunk_size": "library.sidecar:resolve_failure_chunk_size",
}


def __getattr__(name: str) -> Any:
    """Dynamically import exposed helpers and submodules on first access."""

    if name in _LAZY_SUBMODULES:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    if name in _EXPORTS:
        target = _EXPORTS[name]
        module_name, _, attribute = target.partition(":")
        module = import_module(module_name)
        value = getattr(module, attribute) if attribute else module
        globals()[name] = value
        return value
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")


__all__ = sorted((*_LAZY_SUBMODULES, *_EXPORTS))


def __dir__() -> list[str]:
    """Return available attributes including lazily loaded submodules."""

    return sorted({*globals(), *__all__, *_LAZY_SUBMODULES})
