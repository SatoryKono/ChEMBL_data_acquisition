"""Helpers for locating and loading configuration files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

_DEFAULT_CONFIG_NAME = "config.yaml"
_PACKAGE_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = _PACKAGE_ROOT / "config"
DEFAULT_CONFIG_RELATIVE = Path("config") / _DEFAULT_CONFIG_NAME
DEFAULT_CONFIG_PATH = CONFIG_DIR / _DEFAULT_CONFIG_NAME


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
        raise ConfigLoaderError(
            "top-level structure in config file must be a mapping"
        )
    return data, cfg_path


__all__ = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DEFAULT_CONFIG_RELATIVE",
    "ConfigLoaderError",
    "load_yaml_config",
    "resolve_config_path",
]
