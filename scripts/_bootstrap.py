"""Load the project bootstrap helpers when running scripts directly."""

from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path
from typing import Protocol, cast


class _BootstrapModule(Protocol):
    def bootstrap_cli(self, package: str | None, script_file: str | None) -> Path: ...

    def ensure_project_root(self, origin: str | Path | None = None) -> Path: ...


def _project_root_from(script_file: str | None) -> Path | None:
    if not script_file:
        return None
    path = Path(script_file).resolve()
    return path.parent.parent


def _load_bootstrap(script_file: str | None) -> _BootstrapModule:
    module_name = "library.bootstrap"
    module = sys.modules.get(module_name)
    if module is not None:
        return cast(_BootstrapModule, module)

    project_root = _project_root_from(script_file)
    if project_root is not None:
        project_root_str = str(project_root)
        if project_root_str not in sys.path:
            sys.path.insert(0, project_root_str)

    imported = import_module(module_name)
    return cast(_BootstrapModule, imported)


def bootstrap_cli(package: str | None, script_file: str | None) -> None:
    """Ensure the repository root is available for CLI execution."""

    module = _load_bootstrap(script_file)
    module.bootstrap_cli(package, script_file)


def ensure_project_root(origin: str | Path | None = None) -> Path:
    """Expose :func:`library.bootstrap.ensure_project_root` for compatibility."""

    script_hint: str | None
    if origin is None:
        script_hint = None
    else:
        script_hint = str(origin)

    module = _load_bootstrap(script_hint)
    return module.ensure_project_root(origin)


__all__ = ["bootstrap_cli", "ensure_project_root"]
