"""Application configuration using Pydantic models.

This module provides a typed wrapper around ``config.yaml``. Configuration
values are loaded from a YAML file and can be overridden by environment
variables or command line options. The order of precedence is::

    YAML < environment variables < CLI overrides

Environment variables follow the ``CHEMBL_DA__SECTION__KEY`` naming pattern
where sections and keys are joined by double underscores. A number of short
aliases are supported; see ``_ALIAS_MAP`` for the full list.

Additionally, a ``.env`` file located alongside the configuration is loaded
automatically using :mod:`python-dotenv`.
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Sequence
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import yaml
from dotenv import load_dotenv  # type: ignore[import-not-found]
from pydantic import BaseModel, ConfigDict, Field, StrictInt, field_validator
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .log import logger

_EMAIL_RE = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")


def _valid_url(url: str) -> bool:
    """Return ``True`` if *url* contains both scheme and netloc components."""

    parsed = urlparse(url)
    return bool(parsed.scheme and parsed.netloc)


class ConfigError(RuntimeError):
    """Raised when configuration loading fails."""


# ---------------------------------------------------------------------------
# Pydantic models
# ---------------------------------------------------------------------------


class _BaseModel(BaseModel):
    """Base model with common configuration."""

    model_config = ConfigDict(extra="forbid")

    def to_dict(self) -> dict[str, Any]:
        """Return the model as a plain dictionary.

        Uses :meth:`pydantic.BaseModel.model_dump` to obtain a standard
        ``dict`` representation, ensuring compatibility with Pydantic v2.

        Returns
        -------
        dict[str, Any]
            Dictionary representation of the model.


        """

        return self.model_dump()


class _BoolModel(_BaseModel):
    """Model base that parses boolean values from strings."""

    @classmethod
    def _parse_bool(cls, value: Any) -> bool:
        if isinstance(value, bool):
            return value
        val = str(value).strip().lower()
        truthy = {"1", "true", "yes", "on"}
        falsy = {"0", "false", "no", "off"}
        if val in truthy:
            return True
        if val in falsy:
            return False
        raise ValueError(f"Invalid boolean value: {value!r}")


class ApiCfg(_BaseModel):
    """Settings for ChEMBL API access."""

    chembl_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(30, ge=1)
    retries: int = Field(3, ge=0)
    backoff_factor: float = Field(0.5, ge=0)
    rps: int = Field(5, ge=1)
    burst: int = Field(5, ge=1)
    user_agent: str = "chembl-da/0.1 (mailto:info@example.org)"
    api_key: str | None = Field(default=None, repr=False)

    @field_validator("chembl_base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):  # pragma: no cover - defensive
            raise ValueError("invalid URL")
        return v

    @field_validator("user_agent")
    @classmethod
    def _ua(cls, v: str) -> str:
        if not _EMAIL_RE.search(v):
            raise ValueError(
                "api.user_agent must include contact information such as an email"
            )
        return v


class ChemblCfg(_BaseModel):
    cache_ttl: int = Field(3600, ge=1)
    cache_maxsize: int = Field(1024, ge=1)


class OpenAlexCfg(_BaseModel):
    base: str = "https://api.openalex.org"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(30, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)
    mailto: str = "info@example.org"

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    @field_validator("mailto")
    @classmethod
    def _email(cls, v: str) -> str:
        if not v or not _EMAIL_RE.fullmatch(v):
            raise ValueError("openalex.mailto must be a valid email")
        return v


class CrossRefCfg(_BaseModel):
    base: str = "https://api.crossref.org"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(30, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)
    mailto: str = "info@example.org"

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    @field_validator("mailto")
    @classmethod
    def _email(cls, v: str) -> str:
        if not v or not _EMAIL_RE.fullmatch(v):
            raise ValueError("crossref.mailto must be a valid email")
        return v


class UniprotCfg(_BaseModel):
    base: str = "https://rest.uniprot.org"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(30, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)
    delay: float = Field(0.25, ge=0)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class UniprotMappingCfg(_BaseModel):
    """Settings for the UniProt ID mapping service."""

    base: str = "https://rest.uniprot.org/idmapping"
    poll_interval: float = Field(0.5, gt=0)
    timeout: float = Field(300.0, ge=1)
    cache_maxsize: int = Field(128, ge=0)
    cache_ttl: float | None = Field(None, ge=0)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    def __hash__(self) -> int:  # pragma: no cover - simple tuple hash
        """Return a hash based on the configuration fields.

        The model is treated as immutable for caching purposes.
        """

        return hash(
            (
                self.base,
                self.poll_interval,
                self.timeout,
                self.cache_maxsize,
                self.cache_ttl,
            )
        )


class IupharCfg(_BaseModel):
    base: str = "https://www.guidetopharmacology.org/services"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(30, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class PubChemCfg(_BaseModel):
    base: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(60, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(3, ge=1)
    burst: int = Field(5, ge=1)
    delay: float = Field(3.0, ge=0)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class PubMedCfg(_BaseModel):
    base: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(10, ge=1)
    retries: int = Field(2, ge=0)
    encodings: list[str] = Field(
        default_factory=lambda: ["utf-8-sig", "cp1251", "latin1"]
    )

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class SemanticScholarCfg(_BaseModel):
    base: str = "https://api.semanticscholar.org/graph/v1"
    timeout_connect: int = Field(5, ge=1)
    timeout_read: int = Field(10, ge=1)
    retries: int = Field(2, ge=0)
    encodings: list[str] = Field(default_factory=lambda: ["utf-8-sig"])

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class DocTypeCfg(_BaseModel):
    weights: dict[str, int] = Field(
        default_factory=lambda: {"pubmed": 4, "openalex": 3, "scholar": 2}
    )
    thresholds: dict[str, int] = Field(
        default_factory=lambda: {"review": 1, "experimental": 1, "unknown": 2}
    )


class ResourcesCfg(_BaseModel):
    dictionary_dir: Path = Path("dictionary")
    iuphar_target_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_target.csv")
    iuphar_family_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_family.csv")
    uniprot_data_dir: Path = Path("uniprot")
    organism_csv: Path = Path("dictionary/organism.csv")
    status_csv: Path = Path("dictionary/status.csv")
    targets_type_csv: Path = Path("dictionary/targets_type.csv")


class IoCfg(_BoolModel):
    output_dir: Path = Path("data/output")
    cache_dir: Path = Path(".cache")
    csv_sep: str = ","
    csv_encoding: str = "utf-8-sig"
    na_markers: Sequence[str] | None = ("#N/A",)
    csv_index: bool = False
    exist_ok: bool = True

    @field_validator("csv_index", "exist_ok", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class JobsCfg(_BaseModel):
    concurrency: int = Field(8, ge=1)
    chunk_size: int = Field(500, ge=1)


class LogCfg(_BaseModel):
    level: str = "INFO"
    format: str = "[%(asctime)s] %(levelname)s %(name)s: %(message)s"
    datefmt: str = "%Y-%m-%d %H:%M:%S"

    @field_validator("level")
    @classmethod
    def _level(cls, v: str) -> str:
        mapping = getattr(logging, "getLevelNamesMapping", None)
        if callable(mapping):
            valid = {name.upper() for name in mapping()}
        else:  # pragma: no cover - fallback
            valid = {name.upper() for name in logging._nameToLevel}
        if v.upper() not in valid:
            choices = ", ".join(sorted(valid))
            raise ValueError(f"log.level must be one of {choices}, got {v!r}")
        return v


class InitCfg(_BaseModel):
    same_doc: Path = Path("data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx")
    all_doc: Path = Path("data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx")
    output_dir: Path = Path("data/output/ChEMBL/processed")


class BatchCfg(_BoolModel):
    size: int = Field(1000, ge=1)
    pause: float = Field(0.0, ge=0)
    concurrency: int = Field(2, ge=1)
    fail_fast: bool = True
    retry_failed: bool = True

    @field_validator("fail_fast", "retry_failed", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class QualityCfg(_BaseModel):
    sample_rows: int = Field(0, ge=0)
    corr_method: str = "pearson"
    max_unique_preview: int = Field(50, ge=1)
    bin_count: int = Field(20, ge=1)


class MapperCfg(_BoolModel):
    enable_cache: bool = True
    strict_schema: bool = True
    warn_on_cast: bool = True

    @field_validator("enable_cache", "strict_schema", "warn_on_cast", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class RateCfg(_BaseModel):
    global_rps: int = Field(8, ge=1)
    global_burst: int = Field(8, ge=1)
    limiter_cache_maxsize: int = Field(128, ge=1)
    limiter_cache_ttl: int = Field(600, ge=1)


class RetryCfg(_BaseModel):
    max_attempts: int = Field(3, ge=1)
    backoff_factor: float = Field(0.5, ge=0)
    status_forcelist: list[StrictInt] = Field(
        default_factory=lambda: [429, 500, 502, 503, 504]
    )


class ActivityCfg(_BoolModel):
    column: str = "activity_id"
    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = None
    dry_run: bool = False

    @field_validator("dry_run", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class AssayCfg(_BaseModel):
    column: str = "assay_chembl_id"
    chunk_size: int = Field(10, ge=1)
    timeout: float = Field(30.0, ge=0)


class TestitemCfg(_BaseModel):
    column: str = "molecule_chembl_id"
    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)


class DocumentPubmedCfg(_BaseModel):
    column: str = "PMID"
    sleep: float = Field(5.0, ge=0)
    workers: int = Field(1, ge=1)
    batch_size: int = Field(100, ge=1)


class DocumentChemblCfg(_BaseModel):
    column: str = "document_chembl_id"
    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)


class DocumentAllCfg(_BaseModel):
    column: str = "document_chembl_id"
    chunk_size: int = Field(5, ge=1)
    sleep: float = Field(5.0, ge=0)
    workers: int = Field(1, ge=1)
    batch_size: int = Field(50, ge=1)
    timeout: float = Field(30.0, ge=0)


class DocumentCfg(_BaseModel):
    pubmed: DocumentPubmedCfg = Field(default_factory=lambda: DocumentPubmedCfg())
    chembl: DocumentChemblCfg = Field(default_factory=lambda: DocumentChemblCfg())
    all: DocumentAllCfg = Field(default_factory=lambda: DocumentAllCfg())


class TargetUniprotCfg(_BaseModel):
    column: str = "uniprot_id"
    data_dir: Path = Path("dictionary/uniprot")


class TargetChemblCfg(_BaseModel):
    column: str = "chembl_id"
    timeout: float = Field(30.0, ge=0)


class TargetIupharCfg(_BaseModel):
    target_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_target.csv")
    family_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_family.csv")


class TargetAllCfg(_BaseModel):
    data_dir: Path = Path("dictionary/uniprot")
    target_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_target.csv")
    family_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_family.csv")
    timeout: float = Field(30.0, ge=0)
    organism_csv: Path = Path("dictionary/_Target/organism.csv")
    uniprot_column: str = "uniprot_id"
    chembl_out: Path | None = None
    uniprot_out: Path | None = None
    iuphar_out: Path | None = None


class TargetCfg(_BaseModel):
    uniprot: TargetUniprotCfg = Field(default_factory=lambda: TargetUniprotCfg())
    chembl: TargetChemblCfg = Field(default_factory=lambda: TargetChemblCfg())
    iuphar: TargetIupharCfg = Field(default_factory=lambda: TargetIupharCfg())
    all: TargetAllCfg = Field(default_factory=lambda: TargetAllCfg())


class Config(_BaseModel):
    api: ApiCfg = Field(default_factory=lambda: ApiCfg())
    chembl: ChemblCfg = Field(default_factory=lambda: ChemblCfg())
    openalex: OpenAlexCfg = Field(default_factory=lambda: OpenAlexCfg())
    crossref: CrossRefCfg = Field(default_factory=lambda: CrossRefCfg())
    uniprot: UniprotCfg = Field(default_factory=lambda: UniprotCfg())
    uniprot_mapping: UniprotMappingCfg = Field(
        default_factory=lambda: UniprotMappingCfg()
    )
    iuphar: IupharCfg = Field(default_factory=lambda: IupharCfg())
    pubchem: PubChemCfg = Field(default_factory=lambda: PubChemCfg())
    pubmed: PubMedCfg = Field(default_factory=lambda: PubMedCfg())
    semantic_scholar: SemanticScholarCfg = Field(
        default_factory=lambda: SemanticScholarCfg()
    )
    doc_type: DocTypeCfg = Field(default_factory=lambda: DocTypeCfg())
    resources: ResourcesCfg = Field(default_factory=lambda: ResourcesCfg())
    io: IoCfg = Field(default_factory=lambda: IoCfg())
    jobs: JobsCfg = Field(default_factory=lambda: JobsCfg())
    log: LogCfg = Field(default_factory=lambda: LogCfg())
    init: InitCfg = Field(default_factory=lambda: InitCfg())
    batch: BatchCfg = Field(default_factory=lambda: BatchCfg())
    quality: QualityCfg = Field(default_factory=lambda: QualityCfg())
    mapper: MapperCfg = Field(default_factory=lambda: MapperCfg())
    rate: RateCfg = Field(default_factory=lambda: RateCfg())
    retry: RetryCfg = Field(default_factory=lambda: RetryCfg())
    activity: ActivityCfg = Field(default_factory=lambda: ActivityCfg())
    assay: AssayCfg = Field(default_factory=lambda: AssayCfg())
    testitem: TestitemCfg = Field(default_factory=lambda: TestitemCfg())
    document: DocumentCfg = Field(default_factory=lambda: DocumentCfg())
    target: TargetCfg = Field(default_factory=lambda: TargetCfg())


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------


def session_with_retry(api: ApiCfg, retry: RetryCfg) -> Session:
    """Return an HTTP session configured for retries and user agent.

    The returned session retries failed requests for *all* HTTP methods,
    including ``POST``. It also avoids raising exceptions on HTTP error status
    codes, allowing callers to handle responses manually.

    Parameters
    ----------
    api:
        API configuration containing the user agent string.
    retry:
        Retry configuration with maximum attempts and backoff strategy.

    Returns
    -------
    Session
        A ``requests.Session`` instance with retry logic applied.
    """

    session = Session()
    retry_cfg = Retry(
        total=retry.max_attempts,
        backoff_factor=retry.backoff_factor,
        status_forcelist=retry.status_forcelist,
        # ``None`` disables method filtering and retries all HTTP methods.
        # Using ``None`` directly avoids ``Collection`` union type evaluation
        # issues under Python 3.12.
        allowed_methods=None,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry_cfg)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers["User-Agent"] = api.user_agent
    return session


def _set_by_path(data: dict[str, Any], path: list[str], value: Any) -> None:
    cur: dict[str, Any] = data
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[path[-1]] = value


def _apply_env_overrides(data: dict[str, Any]) -> None:
    prefix = "CHEMBL_DA"
    for env_key, env_val in os.environ.items():
        key = env_key.upper()
        if key in _ALIAS_MAP:
            path = _ALIAS_MAP[key]
        elif key.startswith(prefix + "__"):
            path = key[len(prefix) + 2 :].split("__")
        else:
            continue
        parts = [p.lower() for p in path]
        if not _is_valid_path(parts):
            logger.warning(f"Environment variable {key} ignored")
            continue
        _set_by_path(data, parts, env_val)


def _is_valid_path(path: list[str]) -> bool:
    model: type[BaseModel] | None = Config
    for part in path:
        if model is None or part not in model.model_fields:
            return False
        field = model.model_fields[part].annotation
        model = (
            field if isinstance(field, type) and issubclass(field, BaseModel) else None
        )
    return True


def _collect_unknown_keys(
    data: dict[str, Any], model: type[BaseModel], prefix: str = ""
) -> list[str]:
    unknown: list[str] = []
    for key, val in list(data.items()):
        if key not in model.model_fields:
            unknown.append(prefix + key)
            del data[key]
            continue
        field = model.model_fields[key]
        submodel = field.annotation
        if (
            isinstance(val, dict)
            and isinstance(submodel, type)
            and issubclass(submodel, BaseModel)
        ):
            unknown.extend(
                _collect_unknown_keys(val, submodel, prefix=f"{prefix}{key}.")
            )
    return unknown


def load_config(
    path: str | Path = "config.yaml",
    cli_overrides: dict[str, Any] | None = None,
    *,
    strict: bool = False,
) -> Config:
    """Load configuration from *path* applying overrides."""

    load_dotenv(Path(path).resolve().with_name(".env"))

    try:
        with Path(path).open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}
    except FileNotFoundError as exc:  # pragma: no cover - defensive
        raise ConfigError(f"configuration file not found: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(
            f"failed to parse YAML configuration at {path}: {exc}"
        ) from exc

    if not isinstance(data, dict):
        raise TypeError("top-level structure in config file must be a mapping")

    # Guard against accidentally passing the JSON schema instead of a runtime
    # configuration file. The schema contains the ``$defs`` key at the top
    # level which is not present in actual application settings.
    if "$defs" in data:
        raise ConfigError(
            f"{path} appears to be a configuration schema; "
            "provide an application config file such as config.yaml."
        )

    _apply_env_overrides(data)

    if cli_overrides:
        for key, val in cli_overrides.items():
            _set_by_path(data, key.split("."), val)

    unknown = _collect_unknown_keys(data, Config)
    if unknown:
        msg = f"Unknown configuration key(s) in {path}: {', '.join(sorted(unknown))}"
        if strict:
            raise ValueError(msg)
        logger.warning(msg)

    cfg = Config.model_validate(data)

    if not cfg.io.exist_ok:
        for p in (cfg.io.output_dir, cfg.io.cache_dir):
            if not p.exists():
                raise FileNotFoundError(p)

    from .rate_limiter import configure_limiter_cache

    configure_limiter_cache(cfg.rate.limiter_cache_maxsize, cfg.rate.limiter_cache_ttl)
    return cfg


def ensure_dirs(cfg: Config) -> None:
    """Create output and cache directories if required."""

    for path in (cfg.io.output_dir, cfg.io.cache_dir):
        if path.exists():
            if not path.is_dir():
                raise NotADirectoryError(f"{path} is not a directory")
        else:
            if cfg.io.exist_ok:
                path.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(f"{path} does not exist")


def _serialize_paths(data: Any) -> Any:
    if isinstance(data, dict):
        return {k: _serialize_paths(v) for k, v in data.items()}
    if isinstance(data, list):
        return [_serialize_paths(v) for v in data]
    if isinstance(data, Path):
        return str(data)
    return data


def _mask_secrets(data: Any) -> Any:
    secret_tokens = {"key", "token", "secret", "password"}
    if isinstance(data, dict):
        return {
            k: (
                "***"
                if any(t in k.lower() for t in secret_tokens)
                else _mask_secrets(v)
            )
            for k, v in data.items()
        }
    if isinstance(data, list):
        return [_mask_secrets(v) for v in data]
    return data


def print_config(cfg: Config) -> None:
    """Print ``cfg`` as YAML masking secret values."""

    data = _serialize_paths(cfg.to_dict())
    masked = _mask_secrets(data)
    print(yaml.safe_dump(masked, sort_keys=False))


def build_alias_map(
    model: type[BaseModel], prefix: str = "CHEMBL_DA"
) -> dict[str, list[str]]:
    """Generate environment variable aliases for a Pydantic model.

    Parameters
    ----------
    model:
        Root Pydantic model to inspect.
    prefix:
        Alias prefix, defaults to ``"CHEMBL_DA"``.

    Returns
    -------
    dict[str, list[str]]
        Mapping of alias names to configuration paths.
    """

    mapping: dict[str, list[str]] = {}

    def _walk(cls: type[BaseModel], path: list[str]) -> None:
        for name, field in cls.model_fields.items():
            sub_path = path + [name]
            annotation = field.annotation
            if isinstance(annotation, type) and issubclass(annotation, BaseModel):
                _walk(annotation, sub_path)
            else:
                alias = prefix + "_" + "_".join(p.upper() for p in sub_path)
                mapping[alias] = sub_path

    _walk(model, [])
    return mapping


_ALIAS_OVERRIDES: dict[str, list[str]] = {
    "CHEMBL_DA_BASE": ["api", "chembl_base"],
    "CHEMBL_DA_BURST": ["api", "burst"],
    "CHEMBL_DA_CACHE_DIR": ["io", "cache_dir"],
    "CHEMBL_DA_CACHE_MAXSIZE": ["chembl", "cache_maxsize"],
    "CHEMBL_DA_CACHE_TTL": ["chembl", "cache_ttl"],
    "CHEMBL_DA_CHUNK_SIZE": ["jobs", "chunk_size"],
    "CHEMBL_DA_CONCURRENCY": ["jobs", "concurrency"],
    "CHEMBL_DA_DICT_DIR": ["resources", "dictionary_dir"],
    "CHEMBL_DA_GLOBAL_BURST": ["rate", "global_burst"],
    "CHEMBL_DA_GLOBAL_RPS": ["rate", "global_rps"],
    "CHEMBL_DA_IUPHAR_FAMILY_CSV": ["resources", "iuphar_family_csv"],
    "CHEMBL_DA_IUPHAR_TARGET_CSV": ["resources", "iuphar_target_csv"],
    "CHEMBL_DA_LIMITER_CACHE_MAXSIZE": ["rate", "limiter_cache_maxsize"],
    "CHEMBL_DA_LIMITER_CACHE_TTL": ["rate", "limiter_cache_ttl"],
    "CHEMBL_DA_ORGANISM_CSV": ["resources", "organism_csv"],
    "CHEMBL_DA_OUTDIR": ["io", "output_dir"],
    "CHEMBL_DA_RPS": ["api", "rps"],
    "CHEMBL_DA_STATUS_CSV": ["resources", "status_csv"],
    "CHEMBL_DA_TARGETS_TYPE_CSV": ["resources", "targets_type_csv"],
    "CHEMBL_DA_TIMEOUT_CONNECT": ["api", "timeout_connect"],
    "CHEMBL_DA_TIMEOUT_READ": ["api", "timeout_read"],
    "CHEMBL_DA_UNIPROT_DATA_DIR": ["resources", "uniprot_data_dir"],
}

_ALIAS_MAP: dict[str, list[str]] = {
    **build_alias_map(Config),
    **_ALIAS_OVERRIDES,
}


__all__ = [
    "ApiCfg",
    "ChemblCfg",
    "OpenAlexCfg",
    "CrossRefCfg",
    "UniprotCfg",
    "UniprotMappingCfg",
    "IupharCfg",
    "PubChemCfg",
    "PubMedCfg",
    "SemanticScholarCfg",
    "DocTypeCfg",
    "ActivityCfg",
    "AssayCfg",
    "TestitemCfg",
    "DocumentPubmedCfg",
    "DocumentChemblCfg",
    "DocumentAllCfg",
    "DocumentCfg",
    "TargetUniprotCfg",
    "TargetChemblCfg",
    "TargetIupharCfg",
    "TargetAllCfg",
    "TargetCfg",
    "ResourcesCfg",
    "IoCfg",
    "JobsCfg",
    "BatchCfg",
    "QualityCfg",
    "MapperCfg",
    "InitCfg",
    "RateCfg",
    "RetryCfg",
    "session_with_retry",
    "LogCfg",
    "Config",
    "ConfigError",
    "ensure_dirs",
    "load_config",
    "print_config",
    "build_alias_map",
]
