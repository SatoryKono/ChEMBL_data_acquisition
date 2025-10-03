"""Runtime helpers for script execution.

This module ensures that executing a script via ``python scripts/foo.py``
behaves the same as ``python -m scripts.foo`` by pre-pending the repository
root (``library/utils/../..``) to :mod:`sys.path`.
"""

from __future__ import annotations

import sys
from pathlib import Path
from types import ModuleType
from typing import Iterable

_PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _iter_module_paths(module: ModuleType) -> Iterable[Path]:
    """Yield filesystem paths that define ``module``."""

    file_attr = getattr(module, "__file__", None)
    if file_attr:
        yield Path(file_attr).resolve()

    package_paths = getattr(module, "__path__", None)
    if package_paths is not None:
        for package_path in package_paths:
            yield Path(package_path).resolve()


def _is_within_project(path: Path) -> bool:
    """Return ``True`` when ``path`` belongs to the repository tree."""

    try:
        path.relative_to(_PROJECT_ROOT)
    except ValueError:
        return False
    return True


def _should_purge_existing(module: ModuleType | None) -> bool:
    """Check whether ``module`` was loaded from outside the project."""

    if module is None:
        return False
    for module_path in _iter_module_paths(module):
        if _is_within_project(module_path):
            return False
    return True


def _purge_conflicting_modules() -> None:
    """Remove third-party ``library`` modules shadowing the local package."""

    existing = sys.modules.get("library")
    if not _should_purge_existing(existing):
        return

    for name in list(sys.modules):
        if name == "library" or name.startswith("library."):
            del sys.modules[name]


def ensure_project_root() -> None:
    """Insert the project root into :data:`sys.path` and purge conflicts."""

    project_root = str(_PROJECT_ROOT)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
    _purge_conflicting_modules()
