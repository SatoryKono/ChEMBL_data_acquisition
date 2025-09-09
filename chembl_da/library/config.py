"""YAML-based configuration handling for ``chembl_da``.

This module provides a typed configuration system for command-line
utilities. Settings are loaded from ``config.yaml`` and may be overridden
via environment variables.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
from typing import Any, Dict

import yaml

# ---------------------------------------------------------------------------
# Dataclass definitions
# ---------------------------------------------------------------------------


@dataclass
class APIConfig:
    """Settings for external API access."""

    chembl_base: str
    timeout_connect: int
    timeout_read: int
    rps: int
    burst: int


@dataclass
class IOConfig:
    """Directories used for input/output and caching."""

    output_dir: str
    cache_dir: str


@dataclass
class JobsConfig:
    """Parallel execution parameters."""

    concurrency: int
    chunk_size: int


@dataclass
class Config:
    """Top-level configuration container."""

    api: APIConfig
    io: IOConfig
    jobs: JobsConfig

    # Proxy attributes for convenient access -------------------------------
    @property
    def chembl_base(self) -> str:
        return self.api.chembl_base

    @property
    def rps(self) -> int:
        return self.api.rps

    @property
    def burst(self) -> int:
        return self.api.burst

    @property
    def output_dir(self) -> str:
        return self.io.output_dir

    @property
    def cache_dir(self) -> str:
        return self.io.cache_dir

    @property
    def concurrency(self) -> int:
        return self.jobs.concurrency

    @property
    def chunk_size(self) -> int:
        return self.jobs.chunk_size


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

_DEFAULTS: Dict[str, Dict[str, Any]] = {
    "api": {
        "chembl_base": "https://www.ebi.ac.uk/chembl/api/data",
        "timeout_connect": 5,
        "timeout_read": 30,
        "rps": 5,
        "burst": 5,
    },
    "io": {"output_dir": "data/output", "cache_dir": ".cache"},
    "jobs": {"concurrency": 8, "chunk_size": 500},
}

# Mapping of short environment aliases to configuration keys
_ALIAS_MAP: Dict[str, tuple[str, str]] = {
    "CHEMBL_DA_OUTDIR": ("io", "output_dir"),
    "CHEMBL_DA_CACHEDIR": ("io", "cache_dir"),
    "CHEMBL_DA_RPS": ("api", "rps"),
    "CHEMBL_DA_BURST": ("api", "burst"),
    "CHEMBL_DA_TIMEOUT_CONNECT": ("api", "timeout_connect"),
    "CHEMBL_DA_TIMEOUT_READ": ("api", "timeout_read"),
    "CHEMBL_DA_BASE": ("api", "chembl_base"),
    "CHEMBL_DA_CONCURRENCY": ("jobs", "concurrency"),
    "CHEMBL_DA_CHUNK_SIZE": ("jobs", "chunk_size"),
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _deep_update(base: Dict[str, Any], updates: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``updates`` into ``base``."""
    for key, value in updates.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_update(base[key], value)
        else:
            base[key] = value
    return base


def _cast(value: str, current: Any) -> Any:
    """Cast environment ``value`` to the type of ``current``."""
    target = type(current)
    if target is int:
        return int(value)
    if target is bool:
        val = value.lower()
        if val in {"1", "true", "yes"}:
            return True
        if val in {"0", "false", "no"}:
            return False
        raise ValueError(f"invalid boolean value: {value}")
    return value


def _apply_env_hierarchical(cfg: Dict[str, Dict[str, Any]]) -> None:
    """Apply hierarchical ``CHEMBL_DA__SECTION__KEY`` overrides."""
    prefix = "CHEMBL_DA__"
    for env, value in os.environ.items():
        if not env.startswith(prefix):
            continue
        parts = env.split("__", 2)
        if len(parts) != 3:
            continue
        section, key = parts[1].lower(), parts[2].lower()
        if section in cfg and key in cfg[section]:
            cfg[section][key] = _cast(value, cfg[section][key])


def _apply_env_aliases(cfg: Dict[str, Dict[str, Any]]) -> None:
    """Apply short alias environment variable overrides."""
    for env, (section, key) in _ALIAS_MAP.items():
        if env in os.environ:
            cfg[section][key] = _cast(os.environ[env], cfg[section][key])


def _validate(cfg: Config) -> None:
    """Validate the configuration values.

    Raises
    ------
    ValueError
        If any configuration value is invalid.
    """
    if not cfg.api.chembl_base or not cfg.api.chembl_base.startswith(
        ("http://", "https://")
    ):
        raise ValueError("api.chembl_base must start with http:// or https://")
    if cfg.api.timeout_connect <= 0:
        raise ValueError("api.timeout_connect must be positive")
    if cfg.api.timeout_read <= 0:
        raise ValueError("api.timeout_read must be positive")
    if cfg.api.rps <= 0:
        raise ValueError("api.rps must be positive")
    if cfg.api.burst <= 0:
        raise ValueError("api.burst must be positive")
    if cfg.io.output_dir == "":
        raise ValueError("io.output_dir must not be empty")
    if cfg.io.cache_dir == "":
        raise ValueError("io.cache_dir must not be empty")
    if cfg.jobs.concurrency <= 0:
        raise ValueError("jobs.concurrency must be positive")
    if cfg.jobs.chunk_size <= 0:
        raise ValueError("jobs.chunk_size must be positive")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(path: str = "config.yaml") -> Config:
    """Load configuration from ``path`` and environment variables.

    Parameters
    ----------
    path:
        Location of the YAML configuration file. If the file does not exist,
        built-in defaults are used.

    Returns
    -------
    Config
        Parsed configuration object.

    Raises
    ------
    ValueError
        If configuration values are invalid or environment overrides are
        malformed.
    """
    config_dict: Dict[str, Dict[str, Any]] = {
        section: values.copy() for section, values in _DEFAULTS.items()
    }
    config_path = Path(path)
    if config_path.exists():
        data = yaml.safe_load(config_path.read_text()) or {}
        if not isinstance(data, dict):
            raise ValueError(
                "configuration file must contain a mapping at the top level"
            )
        _deep_update(config_dict, data)

    _apply_env_hierarchical(config_dict)
    _apply_env_aliases(config_dict)

    cfg = Config(
        api=APIConfig(**config_dict["api"]),
        io=IOConfig(**config_dict["io"]),
        jobs=JobsConfig(**config_dict["jobs"]),
    )
    _validate(cfg)
    return cfg
