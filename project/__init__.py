"""Compatibility namespace for legacy `project` imports.

Historically the duplicate analysis helpers and other utilities lived under a
`project` package.  The refactor to the consolidated :mod:`library` package
removed that directory which broke older scripts still importing from
``project.*``.  This module keeps those imports working by lazily redirecting
lookups to the new modules.

Use :func:`resolve_legacy_module` to inspect where an attribute will be
forwarded.  New code should import from :mod:`library` directly.
"""

from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import Final

_REDIRECTS: Final[dict[str, str]] = {
    "duplicate_analysis": "project.duplicate_analysis",
}


def resolve_legacy_module(name: str) -> str:
    """Return the fully-qualified module path for *name*.

    Parameters
    ----------
    name:
        Attribute requested on the :mod:`project` package.

    Returns
    -------
    str
        The module that will be imported on demand.

    Raises
    ------
    AttributeError
        If the attribute is not part of the compatibility map.
    """

    try:
        return _REDIRECTS[name]
    except KeyError as exc:  # pragma: no cover - defensive programming
        raise AttributeError(f"module 'project' has no attribute '{name}'") from exc


_definitely_imported: dict[str, ModuleType] = {}


def __getattr__(name: str) -> ModuleType:
    """Dynamically import legacy modules on first access."""

    target = _REDIRECTS.get(name)
    if target is None:
        raise AttributeError(f"module 'project' has no attribute '{name}'")
    module = import_module(target)
    _definitely_imported[name] = module
    return module


def __dir__() -> list[str]:  # pragma: no cover - reflective helper
    return sorted(set(globals()) | set(_REDIRECTS))


__all__ = ["resolve_legacy_module"]
