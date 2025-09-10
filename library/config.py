"""Application configuration utilities."""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict

import yaml

# Prefix for environment variable overrides
ENV_PREFIX = "CHEMBL_"


@dataclass
class APISettings:
    """Settings related to external API access."""

    timeout: float = 30.0


@dataclass
class Config:
    """Top-level configuration container."""

    api: APISettings = field(default_factory=APISettings)


# Mapping of command-line aliases to nested config paths
_ALIASES: dict[str, tuple[str, str]] = {"timeout": ("api", "timeout")}

_config: Config | None = None


def _apply_override(data: Dict[str, Any], path: tuple[str, ...], value: Any) -> None:
    """Set ``value`` at ``path`` within ``data``.

    Parameters
    ----------
    data:
        Mutable mapping representing the config.
    path:
        Sequence of keys describing where to place ``value``.
    value:
        The value to assign.
    """

    cursor = data
    for key in path[:-1]:
        cursor = cursor.setdefault(key, {})  # type: ignore[arg-type]
    cursor[path[-1]] = value


def _build_config(data: Dict[str, Any]) -> Config:
    """Construct :class:`Config` from a nested mapping."""

    api_data = data.get("api", {})
    api = APISettings(**api_data)
    return Config(api=api)


def load_config(path: Path | None = None) -> Config:
    """Load configuration from *path* applying CLI and env overrides.

    Parameters
    ----------
    path:
        Location of the YAML configuration file. Defaults to ``config.yaml`` in
        the project root when ``None``.

    Returns
    -------
    Config
        Parsed configuration with overrides applied.
    """

    cfg_path = path or Path(__file__).resolve().parents[1] / "config.yaml"
    data: Dict[str, Any] = {}
    if cfg_path.exists():
        with cfg_path.open("rt", encoding="utf8") as handle:
            loaded = yaml.safe_load(handle) or {}
            if isinstance(loaded, dict):
                data.update(loaded)

    # Environment variable overrides
    for alias, keys in _ALIASES.items():
        env_name = f"{ENV_PREFIX}{alias.upper()}"
        if (raw := os.getenv(env_name)) is not None:
            _apply_override(data, keys, yaml.safe_load(raw))

    # CLI overrides
    parser = argparse.ArgumentParser(add_help=False)
    for alias in _ALIASES:
        parser.add_argument(f"--{alias}")
    args, _ = parser.parse_known_args()
    for alias, value in vars(args).items():
        if value is not None:
            _apply_override(data, _ALIASES[alias], yaml.safe_load(value))

    return _build_config(data)


def get_config() -> Config:
    """Return the process-wide configuration instance."""

    global _config
    if _config is None:
        _config = load_config()
    return _config


__all__ = ["Config", "load_config", "get_config"]
