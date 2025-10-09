"""Helpers for resolving the installed pipeline version."""

from __future__ import annotations

from functools import lru_cache
from importlib import metadata
from pathlib import Path

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python <3.11
    import tomli as tomllib  # type: ignore[import-not-found, no-redef]

_PACKAGE_NAME = "chembl-data-acquisition"
_DEFAULT_VERSION = "0.0.0"
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
_PYPROJECT_PATH = _PROJECT_ROOT / "pyproject.toml"


@lru_cache(maxsize=1)
def _read_version_from_pyproject() -> str:
    """Return the project version declared in ``pyproject.toml``."""

    if not _PYPROJECT_PATH.exists():
        return _DEFAULT_VERSION
    try:
        with _PYPROJECT_PATH.open("rb") as handle:
            data = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError):  # pragma: no cover - unlikely during tests
        return _DEFAULT_VERSION
    project = data.get("project")
    if isinstance(project, dict):
        version = project.get("version")
        if isinstance(version, str) and version.strip():
            return version
    return _DEFAULT_VERSION


@lru_cache(maxsize=1)
def get_pipeline_version() -> str:
    """Return the installed package version or the value from ``pyproject.toml``."""

    try:
        return metadata.version(_PACKAGE_NAME)
    except metadata.PackageNotFoundError:
        return _read_version_from_pyproject()


__all__ = ["get_pipeline_version"]
