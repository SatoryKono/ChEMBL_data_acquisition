"""Helpers for locating and loading configuration files."""

from __future__ import annotations

import atexit
from contextlib import ExitStack
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

_DEFAULT_CONFIG_NAME = "config.yaml"
_CONFIG_RESOURCE_PACKAGE = "resources.config"
_RESOURCE_STACK = ExitStack()
atexit.register(_RESOURCE_STACK.close)


def _resource_path(*parts: str) -> Path:
    """Return a filesystem path for the packaged config *parts*."""

    traversable = resources.files(_CONFIG_RESOURCE_PACKAGE)
    for part in parts:
        traversable = traversable.joinpath(part)
    return Path(_RESOURCE_STACK.enter_context(resources.as_file(traversable)))


DEFAULT_CONFIG_RELATIVE = Path("resources") / "config" / _DEFAULT_CONFIG_NAME
DEFAULT_CONFIG_PATH = _resource_path(_DEFAULT_CONFIG_NAME)
CONFIG_DIR = DEFAULT_CONFIG_PATH.parent


class ConfigLoaderError(RuntimeError):
    """Raised when configuration parsing fails."""


def resolve_config_path(path: str | Path | None = None) -> Path:
    """Return an absolute path to the configuration file."""

    if path is None:
        return DEFAULT_CONFIG_PATH
    return Path(path).expanduser()


def load_yaml_config(path: str | Path | None = None) -> tuple[dict[str, Any], Path]:
    """Return raw configuration data loaded from YAML."""

    cfg_path = resolve_config_path(path)
    try:
        with cfg_path.open("r", encoding="utf8") as handle:
            data = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:  # pragma: no cover - defensive
        raise ConfigLoaderError(f"configuration file not found: {cfg_path}") from exc
    except yaml.YAMLError as exc:  # pragma: no cover - defensive
        raise ConfigLoaderError(
            f"failed to parse YAML configuration at {cfg_path}: {exc}"
        ) from exc

    if not isinstance(data, dict):
        raise ConfigLoaderError("top-level structure in config file must be a mapping")
    return data, cfg_path


__all__ = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DEFAULT_CONFIG_RELATIVE",
    "ConfigLoaderError",
    "load_yaml_config",
    "resolve_config_path",
]
