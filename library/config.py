"""Configuration utilities for data acquisition scripts.

This module centralises configuration options such as timeouts, rate limits
and output directories. Settings are loaded from ``config.yaml`` when
available, otherwise built-in defaults are used. Values can be overridden
using environment variables prefixed with ``CHEMBL_``. Nested fields are
separated by double underscores, e.g. ``CHEMBL_TIMEOUTS__CONNECT``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
from typing import Any, Dict, Type
from urllib.parse import urlparse

import yaml


class ConfigError(ValueError):
    """Raised when configuration values are invalid."""


@dataclass
class APISettings:
    """Base URLs for external services."""

    chembl_base_url: str = "https://www.ebi.ac.uk/chembl/api/data"
    pubchem_base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    uniprot_base_url: str = "https://rest.uniprot.org"
    eutils_base_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    semanticscholar_base_url: str = "https://api.semanticscholar.org/graph/v1"
    openalex_base_url: str = "https://api.openalex.org"
    crossref_base_url: str = "https://api.crossref.org"
    gtp_base_url: str = "https://www.guidetopharmacology.org"


@dataclass
class TimeoutSettings:
    """Network timeout configuration."""

    connect: float = 10.0
    read: float = 30.0


@dataclass
class RateLimitSettings:
    """HTTP rate limiting and retry settings."""

    max_requests_per_second: float = 5.0
    max_retries: int = 3
    backoff_factor: float = 1.0


@dataclass
class OutputPaths:
    """Output directory configuration."""

    data_dir: Path = Path("data")
    logs_dir: Path = Path("logs")
    tmp_dir: Path = Path("tmp")

    def __post_init__(self) -> None:
        # Allow providing string paths but store as ``Path`` instances.
        self.data_dir = Path(self.data_dir)
        self.logs_dir = Path(self.logs_dir)
        self.tmp_dir = Path(self.tmp_dir)


@dataclass
class Config:
    """Application configuration loaded from ``config.yaml``.

    Parameters
    ----------
    api:
        Base URLs for external services.
    timeouts:
        Network timeout configuration.
    rate_limits:
        HTTP rate limiting and retry configuration.
    output:
        Output directory configuration.
    """

    api: APISettings = field(default_factory=APISettings)
    timeouts: TimeoutSettings = field(default_factory=TimeoutSettings)
    rate_limits: RateLimitSettings = field(default_factory=RateLimitSettings)
    output: OutputPaths = field(default_factory=OutputPaths)

    def __post_init__(self) -> None:
        """Validate configuration values.

        Raises
        ------
        ConfigError
            If any configuration value is invalid.
        """

        # Validate timeouts
        for name, value in vars(self.timeouts).items():
            if value <= 0:
                raise ConfigError(f"timeout '{name}' must be > 0")

        # Validate rate limits
        if self.rate_limits.max_requests_per_second <= 0:
            raise ConfigError("max_requests_per_second must be > 0")
        if self.rate_limits.max_retries < 0:
            raise ConfigError("max_retries must be >= 0")
        if self.rate_limits.backoff_factor < 0:
            raise ConfigError("backoff_factor must be >= 0")

        # Validate URLs
        for name, url in vars(self.api).items():
            parsed = urlparse(url)
            if not parsed.scheme or not parsed.netloc:
                raise ConfigError(f"invalid URL for {name}: {url}")

        # Validate output directories
        for path in [self.output.data_dir, self.output.logs_dir, self.output.tmp_dir]:
            try:
                path.mkdir(parents=True, exist_ok=True)
            except OSError as exc:  # pragma: no cover - unlikely on CI
                raise ConfigError(f"failed to create directory: {path}") from exc
            if not os.access(path, os.W_OK):
                raise ConfigError(f"directory not writable: {path}")


def _apply_env_overrides(data: Dict[str, Any]) -> Dict[str, Any]:
    """Merge environment variable overrides into ``data``.

    Environment variable names must start with ``CHEMBL_``. Nested configuration
    sections are separated by double underscores. For example::

        CHEMBL_TIMEOUTS__CONNECT=5

    Parameters
    ----------
    data:
        Configuration mapping to update in place.

    Returns
    -------
    dict[str, Any]
        The updated mapping.
    """

    prefix = "CHEMBL_"
    for key, value in os.environ.items():
        if not key.startswith(prefix):
            continue
        path = key[len(prefix) :].lower().split("__")
        current = data
        for part in path[:-1]:
            current = current.setdefault(part, {})
        current[path[-1]] = yaml.safe_load(value)
    return data


def _filter_kwargs(data: Dict[str, Any], cls: Type[Any]) -> Dict[str, Any]:
    """Return mapping of keys in ``data`` that match ``cls`` fields."""

    valid_keys = set(getattr(cls, "__dataclass_fields__", {}))
    return {k: v for k, v in data.items() if k in valid_keys}


def load_config(path: str | Path | None = None) -> Config:
    """Return configuration loaded from ``path`` or defaults.

    Parameters
    ----------
    path:
        Optional path to a YAML configuration file. When omitted,
        ``config.yaml`` located at the repository root is used. Missing files
        result in the default configuration being returned.

    Returns
    -------
    Config
        Parsed configuration object.
    """

    cfg_path = Path(path) if path is not None else Path("config.yaml")
    data: Dict[str, Any] = {}
    if cfg_path.exists():
        with cfg_path.open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}

    data = _apply_env_overrides(data)

    return Config(
        api=APISettings(**_filter_kwargs(data.get("api", {}), APISettings)),
        timeouts=TimeoutSettings(
            **_filter_kwargs(data.get("timeouts", {}), TimeoutSettings)
        ),
        rate_limits=RateLimitSettings(
            **_filter_kwargs(data.get("rate_limits", {}), RateLimitSettings)
        ),
        output=OutputPaths(**_filter_kwargs(data.get("output", {}), OutputPaths)),
    )


DEFAULT_CONFIG = load_config()

__all__ = [
    "Config",
    "ConfigError",
    "load_config",
    "DEFAULT_CONFIG",
    "APISettings",
    "TimeoutSettings",
    "RateLimitSettings",
    "OutputPaths",
]
