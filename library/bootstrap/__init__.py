"""Helpers for bootstrapping the project runtime environment."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from pathlib import Path
from types import ModuleType

_DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[2]


def resolve_project_root(origin: str | Path | None = None) -> Path:
    """Return the repository root for ``origin``.

    Parameters
    ----------
    origin:
        Optional path to a module or package. When ``None`` the bootstrap
        module location is used.
    """

    if origin is None:
        start = _DEFAULT_PROJECT_ROOT
    else:
        path = Path(origin).resolve()
        start = path.parent if path.is_file() else path

    for candidate in (start, *start.parents):
        if (candidate / "pyproject.toml").exists():
            return candidate
    return _DEFAULT_PROJECT_ROOT


def _iter_module_paths(module: ModuleType) -> Iterable[Path]:
    """Yield filesystem paths contributing to ``module``."""

    file_attr = getattr(module, "__file__", None)
    if file_attr:
        yield Path(file_attr).resolve()

    package_paths = getattr(module, "__path__", None)
    if package_paths is not None:
        for package_path in package_paths:
            yield Path(package_path).resolve()


def _is_within_project(path: Path, project_root: Path) -> bool:
    """Return ``True`` when ``path`` is inside ``project_root``."""

    try:
        path.relative_to(project_root)
    except ValueError:
        return False
    return True


def _should_purge_existing(module: ModuleType | None, project_root: Path) -> bool:
    """Check whether an existing ``library`` module conflicts with the project."""

    if module is None:
        return False
    return not any(
        _is_within_project(path, project_root) for path in _iter_module_paths(module)
    )


def _purge_conflicting_modules(project_root: Path) -> None:
    """Remove installed ``library`` packages shadowing the repository copy."""

    existing = sys.modules.get("library")
    if not _should_purge_existing(existing, project_root):
        return

    for name in list(sys.modules):
        if name == "library" or name.startswith("library."):
            del sys.modules[name]


def ensure_project_root(origin: str | Path | None = None) -> Path:
    """Insert the project root into :data:`sys.path`.

    Returns
    -------
    pathlib.Path
        The resolved project root.
    """

    project_root = resolve_project_root(origin)
    project_root_str = str(project_root)
    if project_root_str in sys.path:
        sys.path.remove(project_root_str)
    sys.path.insert(0, project_root_str)
    _purge_conflicting_modules(project_root)
    return project_root


def bootstrap_cli(package: str | None, script_file: str | None) -> Path:
    """Bootstrap CLI execution for ``script_file``.

    Parameters
    ----------
    package:
        Value of ``__package__`` for the caller.
    script_file:
        Value of ``__file__`` for the caller.
    """

    return ensure_project_root(script_file)


__all__ = ["bootstrap_cli", "ensure_project_root", "resolve_project_root"]
