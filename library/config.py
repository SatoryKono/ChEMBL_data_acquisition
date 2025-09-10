"""Configuration utilities for data acquisition scripts.

This module centralises configuration options such as timeouts, rate limits
and output directories. Settings are loaded from ``config.yaml`` when
available, otherwise built-in defaults are used.

The :class:`Config` dataclass performs validation in ``__post_init__`` to
ensure supplied values are sensible.  Invalid values raise
:class:`ConfigError`.

"""

from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
from typing import Any, Dict
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

    data_dir: Path | str = Path("data")
    logs_dir: Path | str = Path("logs")
    tmp_dir: Path | str = Path("tmp")

    def __post_init__(self) -> None:
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
            except OSError as exc:
                raise ConfigError(f"failed to create directory: {path}") from exc
            if not os.access(path, os.W_OK):
                raise ConfigError(f"directory not writable: {path}")


def load_config(path: str | Path | None = None) -> Config:
    """Return configuration loaded from ``path`` or defaults.


    Parameters
    ----------
    path:

        Optional path to a YAML configuration file. When omitted, ``config.yaml``
        located at the repository root is used. Missing files result in the
        default configuration being returned.


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

    return Config(
        api=APISettings(**data.get("api", {})),
        timeouts=TimeoutSettings(**data.get("timeouts", {})),
        rate_limits=RateLimitSettings(**data.get("rate_limits", {})),
        output=OutputPaths(**data.get("output", {})),
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

