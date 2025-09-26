"""Helpers for annotating exported tables with pipeline metadata."""

from __future__ import annotations

from datetime import UTC, datetime
from functools import lru_cache
from importlib import metadata
from pathlib import Path
from typing import Final

import pandas as pd

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python <3.11
    import tomli as tomllib  # type: ignore[import-not-found, no-redef]

_PACKAGE_NAME: Final[str] = "chembl-data-acquisition"
_DEFAULT_VERSION: Final[str] = "0.0.0"
_PROJECT_ROOT: Final[Path] = Path(__file__).resolve().parents[1]
_PYPROJECT_PATH: Final[Path] = _PROJECT_ROOT / "pyproject.toml"


@lru_cache(maxsize=1)
def _read_version_from_pyproject() -> str:
    """Return the project version declared in ``pyproject.toml``."""

    if not _PYPROJECT_PATH.exists():
        return _DEFAULT_VERSION
    try:
        with _PYPROJECT_PATH.open("rb") as handle:
            data = tomllib.load(handle)
    except (
        OSError,
        tomllib.TOMLDecodeError,
    ):  # pragma: no cover - unlikely during tests
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


@lru_cache(maxsize=1)
def get_timestamp_utc() -> str:
    """Return an ISO 8601 timestamp representing the pipeline execution time."""

    return datetime.now(UTC).isoformat()


@lru_cache(maxsize=1)
def pipeline_metadata() -> dict[str, str]:
    """Return a mapping with pipeline metadata columns."""

    return {
        "pipeline_version": get_pipeline_version(),
        "timestamp_utc": get_timestamp_utc(),
    }


def add_pipeline_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` with pipeline metadata columns added."""

    if df.empty:
        # Preserve dtypes while ensuring the columns exist for empty frames.
        metadata_values = pipeline_metadata()
        result = df.copy()
        for column, _value in metadata_values.items():
            result[column] = pd.Series(dtype="string")
        return result

    result = df.copy()
    for column, value in pipeline_metadata().items():
        result[column] = value
    return result


__all__ = [
    "add_pipeline_metadata",
    "get_pipeline_version",
    "get_timestamp_utc",
    "pipeline_metadata",
]
