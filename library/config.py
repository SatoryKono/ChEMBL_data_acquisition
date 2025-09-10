"""Typed configuration loader with environment overrides."""

from __future__ import annotations

from dataclasses import dataclass, field
import dataclasses
from pathlib import Path
from typing import Any, Dict, List
import os
import yaml
from urllib.parse import urlparse


# ---------------------------------------------------------------------------
# Dataclass definitions
# ---------------------------------------------------------------------------


@dataclass
class ApiConfig:
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
class IupharConfig:
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
class IoConfig:
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
    pause: int = 0
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
    required_sheets_same_doc: List[str] = field(
        default_factory=lambda: [
            "assay_step5_same_doc",
            "activity_step5_same_doc",
            "document_step5_same_doc",
            "target_step5_same_doc",
            " molecule_step5_same_doc",
            "pairs_same_doc",
        ]
    )
    required_sheets_all_doc: List[str] = field(
        default_factory=lambda: [
            "assay_step5",
            "activity_step5",
            "document_step5",
            "targets_step5",
            "molecule_step5",
            "step5_pairs",
        ]
    )
    output_dir: str = "data/output/ChEMBL/processed"
    same_doc: str = "data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx"
    all_doc: str = "data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx"
    dictionary_dir: str = "dictionary/"
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
    status_forcelist: List[int] = field(
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
    api: ApiConfig = field(default_factory=ApiConfig)
    openalex: OpenAlexConfig = field(default_factory=OpenAlexConfig)
    crossref: CrossrefConfig = field(default_factory=CrossrefConfig)
    uniprot: UniProtConfig = field(default_factory=UniProtConfig)
    iuphar: IupharConfig = field(default_factory=IupharConfig)
    pubchem: PubChemConfig = field(default_factory=PubChemConfig)
    io: IoConfig = field(default_factory=IoConfig)
    jobs: JobsConfig = field(default_factory=JobsConfig)
    batch: BatchConfig = field(default_factory=BatchConfig)
    quality: QualityConfig = field(default_factory=QualityConfig)
    mapper: MapperConfig = field(default_factory=MapperConfig)
    init: InitConfig = field(default_factory=InitConfig)
    rate: RateConfig = field(default_factory=RateConfig)
    retry: RetryConfig = field(default_factory=RetryConfig)
    log: LogConfig = field(default_factory=LogConfig)

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
# Helper utilities
# ---------------------------------------------------------------------------


def _coerce(value: str, current: Any) -> Any:
    if isinstance(current, bool):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    if isinstance(current, int):
        return int(value)
    if isinstance(current, float):
        return float(value)
    return value


def _update_from_dict(obj: Any, data: Dict[str, Any]) -> None:
    for key, val in data.items():
        if not hasattr(obj, key):
            continue
        current = getattr(obj, key)
        if dataclasses.is_dataclass(current):
            if isinstance(val, dict):
                _update_from_dict(current, val)
            continue
        setattr(obj, key, val)


def _set_by_path(cfg: Config, path: List[str], value: str) -> None:
    obj: Any = cfg
    for name in path[:-1]:
        obj = getattr(obj, name, None)
        if obj is None:
            return
    field_name = path[-1]
    if not hasattr(obj, field_name):
        return
    current = getattr(obj, field_name)
    setattr(obj, field_name, _coerce(value, current))


_ALIAS_MAP: Dict[str, List[str]] = {
    "CHEMBL_DA_OUTDIR": ["io", "output_dir"],
    "CHEMBL_DA_RPS": ["api", "rps"],
    "CHEMBL_DA_BURST": ["api", "burst"],
    "CHEMBL_DA_CONCURRENCY": ["jobs", "concurrency"],
    "CHEMBL_DA_CHUNK_SIZE": ["jobs", "chunk_size"],
}


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def _valid_url(url: str) -> bool:
    parsed = urlparse(url)
    return parsed.scheme in {"http", "https"} and bool(parsed.netloc)


def _validate(cfg: Config) -> None:
    if not _valid_url(cfg.api.chembl_base):
        raise ValueError("invalid ChEMBL base URL")
    sections: list[Any] = [
        cfg.api,
        cfg.openalex,
        cfg.crossref,
        cfg.uniprot,
        cfg.iuphar,
        cfg.pubchem,
    ]
    for section in sections:
        base = getattr(section, "base", None)
        if base and not _valid_url(base):
            raise ValueError(f"invalid URL: {base}")
        if section.timeout_connect <= 0 or section.timeout_read <= 0:
            raise ValueError("timeouts must be positive")
        if section.rps <= 0 or section.burst <= 0:
            raise ValueError("rps and burst must be positive")
    if cfg.jobs.concurrency <= 0 or cfg.jobs.chunk_size <= 0:
        raise ValueError("concurrency and chunk_size must be positive")
    if not cfg.io.output_dir.strip() or not cfg.io.cache_dir.strip():
        raise ValueError("directories must not be empty")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(path: str | Path | None = "config.yaml") -> Config:
    cfg = Config()
    if path and Path(path).is_file():
        with Path(path).open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}
        if isinstance(data, dict):
            _update_from_dict(cfg, data)
    prefix = "CHEMBL_DA"
    for env_key, env_val in os.environ.items():
        key = env_key.upper()
        if key in _ALIAS_MAP:
            _set_by_path(cfg, _ALIAS_MAP[key], env_val)
            continue
        if not key.startswith(prefix + "__"):
            continue
        parts = key.split("__")[1:]  # drop prefix
        path_parts = [p.lower() for p in parts]
        _set_by_path(cfg, path_parts, env_val)
    _validate(cfg)
    return cfg


__all__ = ["Config", "load_config"]
