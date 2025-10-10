"""Helpers for resolving import paths for pipeline configuration."""

from __future__ import annotations

import importlib
from typing import Any, cast

from .types import ImportResolutionError

_DEFAULT_SENTINEL = object()


def import_by_path(
    path: str,
    expected_type: type | tuple[type, ...] | object = _DEFAULT_SENTINEL,
) -> Any:
    """Return the object referenced by ``path``.

    Parameters
    ----------
    path:
        Dotted path in the form ``"package.module:attribute"`` or
        ``"package.module.attribute"``.
    expected_type:
        Optional type or tuple of types that the imported object must be an
        instance of. If omitted, the imported object must be callable.
    """

    module_path: str
    attribute: str | None

    if ":" in path:
        module_path, attribute = path.split(":", 1)
    else:
        parts = path.split(".")
        if len(parts) < 2:
            raise ImportResolutionError(f"Invalid dotted path '{path}'")
        module_path = ".".join(parts[:-1])
        attribute = parts[-1]

    try:
        module = importlib.import_module(module_path)
    except ModuleNotFoundError as exc:
        raise ImportResolutionError(f"Cannot import module '{module_path}'") from exc

    obj: Any
    if attribute is None:
        obj = module
    else:
        try:
            obj = getattr(module, attribute)
        except AttributeError as exc:
            raise ImportResolutionError(
                f"Module '{module_path}' has no attribute '{attribute}'"
            ) from exc

    if expected_type is _DEFAULT_SENTINEL:
        if not callable(obj):
            raise ImportResolutionError(
                f"Imported object '{path}' must be callable, got {type(obj).__name__}"
            )
        return obj

    resolved_expected = cast(type | tuple[type, ...], expected_type)

    if not isinstance(obj, resolved_expected):
        expected_repr = _format_expected_type(resolved_expected)
        raise ImportResolutionError(
            f"Imported object '{path}' must be instance of {expected_repr}, "
            f"got {type(obj).__name__}"
        )

    return obj


def _format_expected_type(expected: type | tuple[type, ...]) -> str:
    if isinstance(expected, tuple):
        names = ", ".join(sorted({t.__name__ for t in expected}))
        return f"({names})"
    if isinstance(expected, type):
        return expected.__name__
    return repr(expected)


def resolve_dotted_path(
    path: str,
    expected_type: type | tuple[type, ...] | object = _DEFAULT_SENTINEL,
) -> Any:
    """Backward compatible alias for :func:`import_by_path`.

    Legacy declarative pipeline configurations may invoke
    :func:`resolve_dotted_path` with an ``expected_type`` argument.  The
    helper mirrors :func:`import_by_path` so that both modern and legacy call
    sites behave identically.
    """

    return import_by_path(path, expected_type=expected_type)


__all__ = ["ImportResolutionError", "import_by_path", "resolve_dotted_path"]
