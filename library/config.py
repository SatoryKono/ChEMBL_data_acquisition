"""Configuration loader with environment and CLI overrides.

This module provides a small typed wrapper around ``config.yaml``. Values are
loaded from the YAML file and can be overridden by environment variables or
command line options. The order of precedence is:

``YAML < environment variables < CLI overrides``.

Environment variables follow the ``CHEMBL_DA__SECTION__KEY`` pattern where
sections and keys are joined by double underscores. Short aliases such as
``CHEMBL_DA_RPS`` are also recognised.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
import os
from pathlib import Path
from typing import Any, Dict, List

import yaml
from urllib.parse import urlparse


# ---------------------------------------------------------------------------
# Dataclass definitions
# ---------------------------------------------------------------------------


@dataclass
class ApiCfg:
    """Settings for ChEMBL API access."""

    chembl_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    backoff_factor: float = 0.5
    rps: int = 5
    burst: int = 5
    user_agent: str = "chembl-da/0.1 (+contact: unset)"


@dataclass
class OpenAlexCfg:
    """Settings for the OpenAlex API."""

    base: str = "https://api.openalex.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = ""


@dataclass
class CrossRefCfg:
    """Settings for the CrossRef API."""

    base: str = "https://api.crossref.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = ""


@dataclass
class UniprotCfg:
    """Settings for the UniProt REST API."""

    base: str = "https://rest.uniprot.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5


@dataclass
class IupharCfg:
    """Settings for the IUPHAR API."""

    base: str = "https://www.guidetopharmacology.org/services"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5


@dataclass
class PubChemCfg:
    """Settings for the PubChem PUG REST API."""

    base: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    timeout_connect: int = 5
    timeout_read: int = 60
    retries: int = 3
    rps: int = 3
    burst: int = 5


@dataclass
class IoCfg:
    """Input/output defaults."""

    output_dir: str = "data/output"
    cache_dir: str = ".cache"
    csv_sep: str = ","
    csv_encoding: str = "utf-8-sig"
    csv_index: bool = False
    exist_ok: bool = True


@dataclass
class JobsCfg:
    """Parallel processing settings."""

    concurrency: int = 8
    chunk_size: int = 500


@dataclass
class LogCfg:
    """Logging configuration."""

    level: str = "INFO"
    format: str = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"
    datefmt: str = "%Y-%m-%d %H:%M:%S"


@dataclass
class InitCfg:
    """Workbook paths for ``get_input_initialisation``."""

    same_doc: str = "data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx"
    all_doc: str = "data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx"
    output_dir: str = "data/output/ChEMBL/processed"


@dataclass
class BatchCfg:
    """Batch processing parameters."""

    size: int = 1000
    pause: float = 0
    concurrency: int = 2
    fail_fast: bool = True
    retry_failed: bool = True


@dataclass
class QualityCfg:
    """Options for table quality analysis."""

    sample_rows: int = 0
    corr_method: str = "pearson"
    max_unique_preview: int = 50
    bin_count: int = 20


@dataclass
class MapperCfg:
    """Mapper tool configuration."""

    enable_cache: bool = True
    strict_schema: bool = True
    warn_on_cast: bool = True


@dataclass
class RateCfg:
    """Global rate limiting settings."""

    global_rps: int = 8
    global_burst: int = 8


@dataclass
class RetryCfg:
    """Retry behaviour for HTTP requests."""

    max_attempts: int = 3
    backoff_factor: float = 0.5
    status_forcelist: List[int] = field(
        default_factory=lambda: [429, 500, 502, 503, 504]
    )


@dataclass
class Config:
    """Aggregate project configuration."""

    api: ApiCfg = field(default_factory=ApiCfg)
    openalex: OpenAlexCfg = field(default_factory=OpenAlexCfg)
    crossref: CrossRefCfg = field(default_factory=CrossRefCfg)
    uniprot: UniprotCfg = field(default_factory=UniprotCfg)
    iuphar: IupharCfg = field(default_factory=IupharCfg)
    pubchem: PubChemCfg = field(default_factory=PubChemCfg)
    io: IoCfg = field(default_factory=IoCfg)
    jobs: JobsCfg = field(default_factory=JobsCfg)
    batch: BatchCfg = field(default_factory=BatchCfg)
    quality: QualityCfg = field(default_factory=QualityCfg)
    mapper: MapperCfg = field(default_factory=MapperCfg)
    init: InitCfg = field(default_factory=InitCfg)
    rate: RateCfg = field(default_factory=RateCfg)
    retry: RetryCfg = field(default_factory=RetryCfg)
    log: LogCfg = field(default_factory=LogCfg)

    def to_dict(self) -> Dict[str, Any]:
        """Return the configuration as a plain dictionary."""

        return asdict(self)

    def to_yaml(self) -> str:
        """Serialise the configuration to a YAML string."""

        return yaml.safe_dump(self.to_dict(), sort_keys=False)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _coerce(value: str, current: Any) -> Any:
    """Coerce ``value`` (from env/CLI) to the type of ``current``."""

    if isinstance(current, bool):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    if isinstance(current, int):
        return int(value)
    if isinstance(current, float):
        return float(value)
    return value


def _update_from_dict(
    obj: Any, data: Dict[str, Any], path: List[str] | None = None
) -> None:
    """Recursively update dataclass ``obj`` with ``data`` validating types."""

    path = [] if path is None else path
    for key, val in data.items():
        if not hasattr(obj, key):
            continue
        current = getattr(obj, key)
        if is_dataclass(current):
            if not isinstance(val, dict):
                joined = ".".join(path + [key])
                raise TypeError(f"{joined} must be a mapping")
            _update_from_dict(current, val, path + [key])
            continue
        expected = type(current)
        if not isinstance(val, expected):
            joined = ".".join(path + [key])
            raise TypeError(f"{joined} must be {expected.__name__}, got {val!r}")
        setattr(obj, key, val)


def _set_by_path(cfg: Config, path: List[str], value: Any) -> None:
    """Set ``value`` at ``path`` inside ``cfg`` with type coercion."""

    obj: Any = cfg
    for name in path[:-1]:
        if not hasattr(obj, name):
            raise KeyError(f"unknown config key: {'.'.join(path)}")
        obj = getattr(obj, name)
    field_name = path[-1]
    if not hasattr(obj, field_name):
        raise KeyError(f"unknown config key: {'.'.join(path)}")
    current = getattr(obj, field_name)
    try:
        if isinstance(value, str):
            value = _coerce(value, current)
        elif not isinstance(value, type(current)):
            raise TypeError
    except Exception as exc:  # pragma: no cover - defensive
        joined = ".".join(path)
        raise TypeError(
            f"{joined} must be {type(current).__name__}, got {value!r}"
        ) from exc
    setattr(obj, field_name, value)


_ALIAS_MAP: Dict[str, List[str]] = {
    "CHEMBL_DA_RPS": ["api", "rps"],
    "CHEMBL_DA_BURST": ["api", "burst"],
    "CHEMBL_DA_OUTDIR": ["io", "output_dir"],
    "CHEMBL_DA_CONCURRENCY": ["jobs", "concurrency"],
    "CHEMBL_DA_CHUNK_SIZE": ["jobs", "chunk_size"],
}


def _apply_env_overrides(cfg: Config) -> None:
    """Apply environment variable overrides to ``cfg``."""

    prefix = "CHEMBL_DA"
    for env_key, env_val in os.environ.items():
        key = env_key.upper()
        if key in _ALIAS_MAP:
            _set_by_path(cfg, _ALIAS_MAP[key], env_val)
            continue
        if not key.startswith(prefix + "__"):
            continue
        parts = key.split("__")[1:]
        path_parts = [p.lower() for p in parts]
        _set_by_path(cfg, path_parts, env_val)


def _valid_url(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def _validate(cfg: Config) -> None:
    """Basic sanity checks for configuration values."""

    if not _valid_url(cfg.api.chembl_base):
        raise ValueError("api.chembl_base must be a valid URL")
    if cfg.api.rps <= 0 or cfg.api.burst <= 0:
        raise ValueError("api.rps and api.burst must be positive")
    services: list[tuple[str, Any]] = [
        ("openalex", cfg.openalex),
        ("crossref", cfg.crossref),
        ("uniprot", cfg.uniprot),
        ("iuphar", cfg.iuphar),
        ("pubchem", cfg.pubchem),
    ]
    for name, service in services:
        if not _valid_url(service.base):
            raise ValueError(f"{name}.base must be a valid URL")
        if service.rps <= 0 or service.burst <= 0:
            raise ValueError(f"{name}.rps and {name}.burst must be positive")
    if cfg.jobs.concurrency <= 0 or cfg.jobs.chunk_size <= 0:
        raise ValueError("jobs.concurrency and jobs.chunk_size must be positive")
    if cfg.batch.size <= 0 or cfg.batch.concurrency <= 0:
        raise ValueError("batch.size and batch.concurrency must be positive")
    if cfg.rate.global_rps <= 0 or cfg.rate.global_burst <= 0:
        raise ValueError("rate.global_rps and rate.global_burst must be positive")
    if cfg.retry.max_attempts <= 0 or cfg.retry.backoff_factor < 0:
        raise ValueError(
            "retry.max_attempts must be positive and backoff_factor non-negative"
        )
    if not cfg.io.output_dir.strip() or not cfg.io.cache_dir.strip():
        raise ValueError("io.output_dir and io.cache_dir must not be empty")
    if not cfg.init.same_doc.strip() or not cfg.init.all_doc.strip():
        raise ValueError("init.same_doc and init.all_doc must not be empty")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(
    path: str | Path = "config.yaml", cli_overrides: Dict[str, Any] | None = None
) -> Config:
    """Load configuration from ``path`` applying environment and CLI overrides.

    Parameters
    ----------
    path:
        Location of the YAML configuration file.
    cli_overrides:
        Mapping of ``"section.key"`` paths to values coming from the command
        line.

    Returns
    -------
    Config
        Fully populated configuration object.
    """

    cfg = Config()
    if path and Path(path).is_file():
        with Path(path).open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}
        if not isinstance(data, dict):
            raise TypeError("top-level structure in config file must be a mapping")
        _update_from_dict(cfg, data)

    _apply_env_overrides(cfg)

    if cli_overrides:
        for key, val in cli_overrides.items():
            _set_by_path(cfg, key.split("."), val)

    _validate(cfg)
    return cfg


__all__ = [
    "ApiCfg",
    "OpenAlexCfg",
    "CrossRefCfg",
    "UniprotCfg",
    "IupharCfg",
    "PubChemCfg",
    "IoCfg",
    "JobsCfg",
    "BatchCfg",
    "QualityCfg",
    "MapperCfg",
    "InitCfg",
    "RateCfg",
    "RetryCfg",
    "LogCfg",
    "Config",
    "load_config",
]
