"""Configuration utilities for data acquisition scripts.

This module centralises configuration handling for the project.  Settings are
loaded from a YAML file and validated using dataclasses.  Values in the YAML can
be overridden by environment variables such as ``CHEMBL_API_BASE_URL`` or
``TIMEOUT_CONNECT``.  Short aliases may also be used in the YAML for brevity,
e.g. ``api.chembl.url`` which maps to ``APISettings.chembl_base_url``.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple, Type
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
        # Accept ``str`` inputs but store as :class:`Path`.
        self.data_dir = Path(self.data_dir)
        self.logs_dir = Path(self.logs_dir)
        self.tmp_dir = Path(self.tmp_dir)


@dataclass
class Config:
    """Validated application configuration.

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

    def __post_init__(self) -> None:  # noqa: D401 - short description inherited
        """Validate configuration values."""

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
            except OSError as exc:  # pragma: no cover - filesystem errors rare
                raise ConfigError(f"failed to create directory: {path}") from exc
            if not os.access(path, os.W_OK):
                raise ConfigError(f"directory not writable: {path}")


# ---------------------------------------------------------------------------
# Aliases and environment variable overrides
# ---------------------------------------------------------------------------

# Mapping of dotted alias paths in the YAML file to canonical configuration
# keys used by the dataclasses.
_ALIAS_MAP: Dict[Tuple[str, ...], Tuple[str, ...]] = {
    ("api", "chembl", "url"): ("api", "chembl_base_url"),
    ("api", "pubchem", "url"): ("api", "pubchem_base_url"),
    ("api", "uniprot", "url"): ("api", "uniprot_base_url"),
    ("api", "eutils", "url"): ("api", "eutils_base_url"),
    ("api", "semanticscholar", "url"): ("api", "semanticscholar_base_url"),
    ("api", "openalex", "url"): ("api", "openalex_base_url"),
    ("api", "crossref", "url"): ("api", "crossref_base_url"),
    ("api", "gtp", "url"): ("api", "gtp_base_url"),
}


# Environment variables overriding configuration values.  Each entry maps the
# variable name to a tuple of (config path, type-caster).
_ENV_VAR_MAP: Dict[str, Tuple[Tuple[str, ...], Type[Any]]] = {
    "CHEMBL_API_BASE_URL": (("api", "chembl_base_url"), str),
    "PUBCHEM_API_BASE_URL": (("api", "pubchem_base_url"), str),
    "UNIPROT_API_BASE_URL": (("api", "uniprot_base_url"), str),
    "EUTILS_API_BASE_URL": (("api", "eutils_base_url"), str),
    "SEMANTICSCHOLAR_API_BASE_URL": (("api", "semanticscholar_base_url"), str),
    "OPENALEX_API_BASE_URL": (("api", "openalex_base_url"), str),
    "CROSSREF_API_BASE_URL": (("api", "crossref_base_url"), str),
    "GTP_API_BASE_URL": (("api", "gtp_base_url"), str),
    "TIMEOUT_CONNECT": (("timeouts", "connect"), float),
    "TIMEOUT_READ": (("timeouts", "read"), float),
    "RATE_LIMIT_MAX_REQUESTS_PER_SECOND": (
        ("rate_limits", "max_requests_per_second"),
        float,
    ),
    "RATE_LIMIT_MAX_RETRIES": (("rate_limits", "max_retries"), int),
    "RATE_LIMIT_BACKOFF_FACTOR": (("rate_limits", "backoff_factor"), float),
    "OUTPUT_DATA_DIR": (("output", "data_dir"), str),
    "OUTPUT_LOGS_DIR": (("output", "logs_dir"), str),
    "OUTPUT_TMP_DIR": (("output", "tmp_dir"), str),
}


def _deep_get(data: Dict[str, Any], path: Tuple[str, ...]) -> Any | None:
    """Retrieve a value from ``data`` following ``path``."""

    cur: Any = data
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return None
        cur = cur[key]
    return cur


def _deep_set(data: Dict[str, Any], path: Tuple[str, ...], value: Any) -> None:
    """Set ``value`` in ``data`` at ``path`` creating nested dictionaries."""

    cur = data
    for key in path[:-1]:
        cur = cur.setdefault(key, {})
    cur[path[-1]] = value


def _apply_aliases(data: Dict[str, Any]) -> None:
    """Expand short-form aliases within ``data``."""

    for alias, target in _ALIAS_MAP.items():
        value = _deep_get(data, alias)
        if value is not None:
            _deep_set(data, target, value)


def _override_with_env(data: Dict[str, Any], env: Mapping[str, str]) -> None:
    """Override ``data`` with values from ``env`` according to ``_ENV_VAR_MAP``."""

    for var, (path, caster) in _ENV_VAR_MAP.items():
        if var in env:
            raw = env[var]
            try:
                value = caster(raw)
            except (TypeError, ValueError):  # pragma: no cover - defensive
                value = raw
            _deep_set(data, path, value)


def _filter_section(data: Dict[str, Any], cls: Type[Any]) -> Dict[str, Any]:
    """Return ``data`` with keys not present on ``cls`` removed."""

    valid = set(cls.__dataclass_fields__.keys())
    return {k: v for k, v in data.items() if k in valid}


def load_config(path: str | Path | None = None) -> Config:
    """Load configuration from ``path`` and apply overrides.

    Parameters
    ----------
    path:
        Optional path to a YAML configuration file. When omitted, ``config.yaml``
        in the repository root is used. Missing files yield the default
        configuration.

    Returns
    -------
    Config
        Parsed configuration object with environment variable overrides
        applied.
    """

    cfg_path = Path(path) if path is not None else Path("config.yaml")
    data: Dict[str, Any] = {}
    if cfg_path.exists():
        with cfg_path.open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}

    _apply_aliases(data)
    _override_with_env(data, os.environ)

    return Config(
        api=APISettings(**_filter_section(data.get("api", {}), APISettings)),
        timeouts=TimeoutSettings(
            **_filter_section(data.get("timeouts", {}), TimeoutSettings)
        ),
        rate_limits=RateLimitSettings(
            **_filter_section(data.get("rate_limits", {}), RateLimitSettings)
        ),
        output=OutputPaths(**_filter_section(data.get("output", {}), OutputPaths)),
    )


DEFAULT_CONFIG = load_config()


__all__ = [
    "Config",
    "ConfigError",
    "APISettings",
    "TimeoutSettings",
    "RateLimitSettings",
    "OutputPaths",
    "load_config",
    "DEFAULT_CONFIG",
]
