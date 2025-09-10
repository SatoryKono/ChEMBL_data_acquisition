"""Typed configuration loader for ChEMBL data acquisition.

This module provides a ``load_config`` function that reads settings from a
YAML file and environment variables.  Values are exposed via dataclasses for
static type checking.  Environment variables take precedence over the YAML
file which in turn overrides the built-in defaults.  If the configuration file
is missing the defaults are used.

Environment variable hierarchy
==============================
Variables use the ``CHEMBL_DA__SECTION__KEY`` form where section and key are
case insensitive::

    export CHEMBL_DA__API__RPS=2

In addition to the hierarchical form a number of short aliases are provided:
``CHEMBL_DA_OUTDIR``, ``CHEMBL_DA_RPS``, ``CHEMBL_DA_BURST``,
``CHEMBL_DA_CONCURRENCY`` and ``CHEMBL_DA_CHUNK_SIZE``.

Examples
========
>>> from chembl_da.library.config import load_config
>>> cfg = load_config("config.yaml")
>>> cfg.api.chembl_base
'https://www.ebi.ac.uk/chembl/api/data'
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Mapping
import os
import yaml
from urllib.parse import urlparse

# ---------------------------------------------------------------------------
# Dataclass definitions mirroring ``config.yaml``


@dataclass
class APIConfig:
    chembl_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    backoff_factor: float = 0.5
    rps: int = 5
    burst: int = 5
    user_agent: str = "chembl-da/0.1 (+contact: unset)"


@dataclass
class OpenAlexConfig:
    base: str = "https://api.openalex.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = ""


@dataclass
class CrossrefConfig:
    base: str = "https://api.crossref.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = ""


@dataclass
class UniProtConfig:
    base: str = "https://rest.uniprot.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5


@dataclass
class IUPHARConfig:
    base: str = "https://www.guidetopharmacology.org/services"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5


@dataclass
class PubChemConfig:
    base: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    timeout_connect: int = 5
    timeout_read: int = 60
    retries: int = 3
    rps: int = 3
    burst: int = 5


@dataclass
class IOConfig:
    output_dir: str = "data/output"
    cache_dir: str = ".cache"
    csv_sep: str = ","
    csv_encoding: str = "utf-8-sig"
    csv_index: bool = False
    xlsx_engine: str = "openpyxl"
    xlsx_min_version: str = "3.1.0"
    exist_ok: bool = True
    compression: str = ""


@dataclass
class JobsConfig:
    concurrency: int = 8
    chunk_size: int = 500


@dataclass
class BatchConfig:
    size: int = 1000
    pause: float = 0.0
    concurrency: int = 2
    fail_fast: bool = True
    retry_failed: bool = True


@dataclass
class QualityConfig:
    sample_rows: int = 0
    corr_method: str = "pearson"
    max_unique_preview: int = 50
    bin_count: int = 20


@dataclass
class MapperConfig:
    enable_cache: bool = True
    strict_schema: bool = True
    warn_on_cast: bool = True


@dataclass
class InitConfig:
    required_sheets: list[str] = field(
        default_factory=lambda: [
            "assay_step5_same_doc",
            "activity_step5_same_doc",
            "documents_step5_same_doc",
            "targets_step5_same_doc",
            "testitem_step5_same_doc",
        ]
    )
    output_dir: str = "data/input/ChEMBL/processed"
    save_pairs_sheets: bool = True
    fail_fast: bool = True


@dataclass
class RateConfig:
    global_rps: int = 8
    global_burst: int = 8


@dataclass
class RetryConfig:
    max_attempts: int = 3
    backoff_factor: float = 0.5
    status_forcelist: list[int] = field(
        default_factory=lambda: [429, 500, 502, 503, 504]
    )


@dataclass
class LogConfig:
    level: str = "INFO"
    format: str = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"
    datefmt: str = "%Y-%m-%d %H:%M:%S"
    color: bool = True


@dataclass
class Config:
    api: APIConfig = field(default_factory=APIConfig)
    openalex: OpenAlexConfig = field(default_factory=OpenAlexConfig)
    crossref: CrossrefConfig = field(default_factory=CrossrefConfig)
    uniprot: UniProtConfig = field(default_factory=UniProtConfig)
    iuphar: IUPHARConfig = field(default_factory=IUPHARConfig)
    pubchem: PubChemConfig = field(default_factory=PubChemConfig)
    io: IOConfig = field(default_factory=IOConfig)
    jobs: JobsConfig = field(default_factory=JobsConfig)
    batch: BatchConfig = field(default_factory=BatchConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
    mapper: MapperConfig = field(default_factory=MapperConfig)
    init: InitConfig = field(default_factory=InitConfig)
    rate: RateConfig = field(default_factory=RateConfig)
    retry: RetryConfig = field(default_factory=RetryConfig)
    log: LogConfig = field(default_factory=LogConfig)

    # Short proxies ---------------------------------------------------------
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
    def concurrency(self) -> int:
        return self.jobs.concurrency

    @property
    def chunk_size(self) -> int:
        return self.jobs.chunk_size


# ---------------------------------------------------------------------------
# Utilities


_ALIAS_MAP: Mapping[str, str] = {
    "CHEMBL_DA_OUTDIR": "io.output_dir",
    "CHEMBL_DA_RPS": "api.rps",
    "CHEMBL_DA_BURST": "api.burst",
    "CHEMBL_DA_CONCURRENCY": "jobs.concurrency",
    "CHEMBL_DA_CHUNK_SIZE": "jobs.chunk_size",
}

_CONFIG: Config | None = None


def _deep_update(base: Dict[str, Any], override: Mapping[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``override`` into ``base``."""
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, Mapping):
            base[key] = _deep_update(dict(base[key]), value)
        else:
            base[key] = value
    return base


def _cast(value: str, target: Any) -> Any:
    """Cast ``value`` to the type of ``target``."""
    if not isinstance(value, str):
        return value
    if isinstance(target, bool):
        return value.lower() in {"1", "true", "yes", "on"}
    if isinstance(target, int):
        return int(value)
    if isinstance(target, float):
        return float(value)
    return value


def _apply_env(data: Dict[str, Any]) -> None:
    """Override ``data`` in-place with environment variables."""
    # Aliases first
    for env, path in _ALIAS_MAP.items():
        if env in os.environ:
            sections = path.split(".")
            ref = data
            for sec in sections[:-1]:
                ref = ref.setdefault(sec, {})
            ref[sections[-1]] = os.environ[env]

    prefix = "CHEMBL_DA__"
    for env, value in os.environ.items():
        if not env.upper().startswith(prefix):
            continue
        parts = env[len(prefix) :].lower().split("__")
        ref = data
        for part in parts[:-1]:
            ref = ref.setdefault(part, {})
        ref[parts[-1]] = value


def _validate(cfg: Config) -> None:
    """Perform basic validation of ``cfg``."""

    def _check_url(url: str) -> None:
        parsed = urlparse(url)
        if not (parsed.scheme and parsed.netloc):
            raise ValueError(f"invalid URL: {url}")

    for section in [
        cfg.api,
        cfg.openalex,
        cfg.crossref,
        cfg.uniprot,
        cfg.iuphar,
        cfg.pubchem,
    ]:
        base = getattr(section, "chembl_base", None) or getattr(section, "base", "")
        _check_url(base)
        if section.timeout_connect <= 0 or section.timeout_read <= 0:
            raise ValueError("timeouts must be positive")
        if section.rps <= 0 or section.burst <= 0:
            raise ValueError("rps and burst must be positive")

    if not cfg.io.output_dir or not cfg.io.cache_dir:
        raise ValueError("I/O directories must not be empty")
    if cfg.jobs.concurrency <= 0 or cfg.jobs.chunk_size <= 0:
        raise ValueError("jobs.concurrency and jobs.chunk_size must be positive")
    if cfg.rate.global_rps <= 0 or cfg.rate.global_burst <= 0:
        raise ValueError("rate limits must be positive")
    if cfg.retry.max_attempts <= 0 or cfg.retry.backoff_factor <= 0:
        raise ValueError("retry settings must be positive")


def load_config(path: str | Path = "config.yaml") -> Config:
    """Load configuration from *path* and environment variables.

    The result is cached so subsequent calls return the same instance unless a
    different ``path`` is supplied.
    """
    global _CONFIG
    if _CONFIG is not None and Path(path) == getattr(_CONFIG, "_source", Path(path)):
        return _CONFIG

    data: Dict[str, Any] = (
        yaml.safe_load(Path(path).read_text()) if Path(path).exists() else {}
    )
    defaults = yaml.safe_load(
        yaml.dump(
            {
                "api": APIConfig().__dict__,
                "openalex": OpenAlexConfig().__dict__,
                "crossref": CrossrefConfig().__dict__,
                "uniprot": UniProtConfig().__dict__,
                "iuphar": IUPHARConfig().__dict__,
                "pubchem": PubChemConfig().__dict__,
                "io": IOConfig().__dict__,
                "jobs": JobsConfig().__dict__,
                "batch": BatchConfig().__dict__,
                "quality": QualityConfig().__dict__,
                "mapper": MapperConfig().__dict__,
                "init": InitConfig().__dict__,
                "rate": RateConfig().__dict__,
                "retry": RetryConfig().__dict__,
                "log": LogConfig().__dict__,
            }
        )
    )
    merged = _deep_update(defaults, data)
    _apply_env(merged)

    def _build(section_cls, values: Mapping[str, Any]):
        inst = section_cls(
            **{k: _cast(v, getattr(section_cls(), k)) for k, v in values.items()}
        )
        return inst

    cfg = Config(
        api=_build(APIConfig, merged.get("api", {})),
        openalex=_build(OpenAlexConfig, merged.get("openalex", {})),
        crossref=_build(CrossrefConfig, merged.get("crossref", {})),
        uniprot=_build(UniProtConfig, merged.get("uniprot", {})),
        iuphar=_build(IUPHARConfig, merged.get("iuphar", {})),
        pubchem=_build(PubChemConfig, merged.get("pubchem", {})),
        io=_build(IOConfig, merged.get("io", {})),
        jobs=_build(JobsConfig, merged.get("jobs", {})),
        batch=_build(BatchConfig, merged.get("batch", {})),
        quality=_build(QualityConfig, merged.get("quality", {})),
        mapper=_build(MapperConfig, merged.get("mapper", {})),
        init=_build(InitConfig, merged.get("init", {})),
        rate=_build(RateConfig, merged.get("rate", {})),
        retry=_build(RetryConfig, merged.get("retry", {})),
        log=_build(LogConfig, merged.get("log", {})),
    )
    _validate(cfg)
    cfg._source = Path(path)  # type: ignore[attr-defined]
    _CONFIG = cfg
    return cfg
