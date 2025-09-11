"""Configuration loader with environment and CLI overrides.

This module provides a small typed wrapper around ``config.yaml``. Values are
loaded from the YAML file and can be overridden by environment variables or
command line options. The order of precedence is:

``YAML < environment variables < CLI overrides``.

Environment variables follow the ``CHEMBL_DA__SECTION__KEY`` pattern where
sections and keys are joined by double underscores. Short aliases such as
``CHEMBL_DA_RPS`` are also recognised. Additional convenience aliases include:

* ``CHEMBL_DA_BASE`` → ``api.chembl_base``
* ``CHEMBL_DA_TIMEOUT_CONNECT`` → ``api.timeout_connect``
* ``CHEMBL_DA_TIMEOUT_READ`` → ``api.timeout_read``

Failure modes
-------------
``load_config`` raises :class:`ConfigError` when the configuration file is
missing or cannot be parsed. Type and value mismatches are reported via the
appropriate built-in exceptions during validation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
import logging
import os
import re
from pathlib import Path
from typing import Any, Dict, List
from urllib.parse import urlparse

import yaml
import jsonschema
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


logger = logging.getLogger(__name__)


_EMAIL_RE = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")


def _valid_url(url: str) -> bool:
    """Check that a URL contains both a scheme and network location.

    Parameters
    ----------
    url:
        URL string to validate.

    Returns
    -------
    bool
        ``True`` when ``url`` contains scheme and netloc components.
    """

    parsed = urlparse(url)
    return bool(parsed.scheme and parsed.netloc)


class ConfigError(RuntimeError):
    """Raised when configuration loading fails."""


# ---------------------------------------------------------------------------
# Dataclass definitions
# ---------------------------------------------------------------------------


@dataclass
class ApiCfg:
    """Settings for ChEMBL API access.

    The ``user_agent`` must include contact information such as an e-mail
    address, e.g. ``"chembl-da/0.1 (mailto:info@example.org)"``.
    """

    chembl_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    backoff_factor: float = 0.5
    rps: int = 5
    burst: int = 5
    user_agent: str = "chembl-da/0.1 (mailto:info@example.org)"


@dataclass
class OpenAlexCfg:
    """Settings for the OpenAlex API.

    ``mailto`` is required by the service and must be a valid e-mail address.
    """

    base: str = "https://api.openalex.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = "info@example.org"


@dataclass
class CrossRefCfg:
    """Settings for the CrossRef API.

    ``mailto`` must be supplied and contain a valid e-mail address.
    """

    base: str = "https://api.crossref.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    mailto: str = "info@example.org"


@dataclass
class UniprotCfg:
    """Settings for the UniProt REST API."""

    base: str = "https://rest.uniprot.org"
    timeout_connect: int = 5
    timeout_read: int = 30
    retries: int = 3
    rps: int = 4
    burst: int = 5
    delay: float = 0.25


@dataclass
class UniprotMappingCfg:
    """Settings for the UniProt ID Mapping API."""

    base: str = "https://rest.uniprot.org/idmapping"
    poll_interval: float = 0.5
    timeout: float = 300.0


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
    delay: float = 3.0


@dataclass
class PubMedCfg:
    """Settings for the PubMed API and related I/O."""

    base: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    timeout_connect: int = 5
    timeout_read: int = 10
    retries: int = 2
    encodings: List[str] = field(
        default_factory=lambda: ["utf-8-sig", "cp1251", "latin1"]
    )


@dataclass
class SemanticScholarCfg:
    """Settings for the Semantic Scholar API."""

    base: str = "https://api.semanticscholar.org/graph/v1"
    timeout_connect: int = 5
    timeout_read: int = 10
    retries: int = 2
    encodings: List[str] = field(default_factory=lambda: ["utf-8-sig"])


@dataclass
class DocTypeCfg:
    """Settings for document type classification."""

    weights: Dict[str, int] = field(
        default_factory=lambda: {"pubmed": 4, "openalex": 3, "scholar": 2}
    )
    thresholds: Dict[str, int] = field(
        default_factory=lambda: {"review": 1, "experimental": 1, "unknown": 2}
    )


@dataclass
class ResourcesCfg:
    """Paths to static resource files used by the application."""

    dictionary_dir: Path = Path("dictionary")
    iuphar_target_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_target.csv")
    iuphar_family_csv: Path = Path("dictionary/_IUPHAR/_IUPHAR_family.csv")
    uniprot_data_dir: Path = Path("uniprot")
    organism_csv: Path = Path("dictionary/organism.csv")
    status_csv: Path = Path("dictionary/status.csv")
    targets_type_csv: Path = Path("dictionary/targets_type.csv")


@dataclass
class IoCfg:
    """Input/output defaults."""

    output_dir: Path = Path("data/output")
    cache_dir: Path = Path(".cache")
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

    same_doc: Path = Path("data/input/ChEMBL/ChEMBL_same_document_20_05.xlsx")
    all_doc: Path = Path("data/input/ChEMBL/ChEMBL_all_10_05_step5.xlsx")
    output_dir: Path = Path("data/output/ChEMBL/processed")


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


def session_with_retry(api: ApiCfg, retry: RetryCfg) -> Session:
    """Return an HTTP session configured for retries and user agent.

    Parameters
    ----------
    api:
        Global API settings providing the ``User-Agent`` header.
    retry:
        Retry configuration specifying backoff behaviour.

    Returns
    -------
    Session
        Configured :class:`requests.Session` instance.
    """

    retry_cfg = Retry(
        total=retry.max_attempts,
        backoff_factor=retry.backoff_factor,
        status_forcelist=retry.status_forcelist,
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry_cfg)
    session = Session()
    session.headers["User-Agent"] = api.user_agent
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


@dataclass
class Config:
    """Aggregate project configuration."""

    api: ApiCfg = field(default_factory=ApiCfg)
    openalex: OpenAlexCfg = field(default_factory=OpenAlexCfg)
    crossref: CrossRefCfg = field(default_factory=CrossRefCfg)
    uniprot: UniprotCfg = field(default_factory=UniprotCfg)
    uniprot_mapping: UniprotMappingCfg = field(default_factory=UniprotMappingCfg)
    iuphar: IupharCfg = field(default_factory=IupharCfg)
    pubchem: PubChemCfg = field(default_factory=PubChemCfg)

    pubmed: PubMedCfg = field(default_factory=PubMedCfg)
    semantic_scholar: SemanticScholarCfg = field(default_factory=SemanticScholarCfg)

    doc_type: DocTypeCfg = field(default_factory=DocTypeCfg)

    resources: ResourcesCfg = field(default_factory=ResourcesCfg)

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
    """Coerce ``value`` (from env/CLI) to the type of ``current``.

    When ``current`` is :class:`bool`, only a restricted set of string
    representations is accepted. Truthy values are ``{"1", "true", "yes",
    "on"}`` and falsy values are ``{"0", "false", "no", "off"}``.
    Any other value raises :class:`ValueError` to avoid silently interpreting
    unexpected input.
    """

    if isinstance(current, bool):
        val = value.strip().lower()
        truthy = {"1", "true", "yes", "on"}
        falsy = {"0", "false", "no", "off"}
        if val in truthy:
            return True
        if val in falsy:
            return False
        raise ValueError(f"Invalid boolean value: {value!r}")
    if isinstance(current, int):
        return int(value)
    if isinstance(current, float):
        return float(value)
    if isinstance(current, Path):
        return Path(value)
    return value


def _update_from_dict(
    obj: Any,
    data: Dict[str, Any],
    path: List[str] | None = None,
    *,
    unknown_keys: List[str] | None = None,
) -> None:
    """Recursively update dataclass ``obj`` with ``data`` validating types.

    Parameters
    ----------
    obj:
        Dataclass instance to update.
    data:
        Mapping of field names to values.
    path:
        Current traversal path used for error reporting.
    unknown_keys:
        Collects dotted-path names of keys not present on ``obj``.
    """

    path = [] if path is None else path
    if unknown_keys is None:
        unknown_keys = []
    for key, val in data.items():
        if not hasattr(obj, key):
            unknown_keys.append(".".join(path + [key]))
            continue
        current = getattr(obj, key)
        if is_dataclass(current):
            if not isinstance(val, dict):
                joined = ".".join(path + [key])
                raise TypeError(f"{joined} must be a mapping")
            _update_from_dict(current, val, path + [key], unknown_keys=unknown_keys)
            continue
        if isinstance(val, str):
            try:
                val = _coerce(val, current)
            except Exception as exc:  # pragma: no cover - defensive
                joined = ".".join(path + [key])
                raise TypeError(
                    f"{joined} must be {type(current).__name__}, got {val!r}"
                ) from exc
        if not isinstance(val, type(current)):
            joined = ".".join(path + [key])
            raise TypeError(f"{joined} must be {type(current).__name__}, got {val!r}")
        setattr(obj, key, val)


def _set_by_path(cfg: Config, path: List[str], value: Any) -> None:
    """Set ``value`` at ``path`` inside ``cfg`` with type coercion."""

    obj: Any = cfg
    for name in path[:-1]:
        if isinstance(obj, dict):
            if name not in obj:
                raise KeyError(f"unknown config key: {'.'.join(path)}")
            obj = obj[name]
            continue
        if not hasattr(obj, name):
            raise KeyError(f"unknown config key: {'.'.join(path)}")
        obj = getattr(obj, name)
    field_name = path[-1]
    if isinstance(obj, dict):
        if field_name not in obj:
            raise KeyError(f"unknown config key: {'.'.join(path)}")
        current = obj[field_name]
    else:
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
        if isinstance(exc, ValueError):
            raise ValueError(f"{joined}: {exc}") from exc
        raise TypeError(
            f"{joined} must be {type(current).__name__}, got {value!r}"
        ) from exc
    if isinstance(obj, dict):
        obj[field_name] = value
    else:
        setattr(obj, field_name, value)


_ALIAS_MAP: Dict[str, List[str]] = {
    "CHEMBL_DA_RPS": ["api", "rps"],
    "CHEMBL_DA_BURST": ["api", "burst"],
    "CHEMBL_DA_BASE": ["api", "chembl_base"],
    "CHEMBL_DA_TIMEOUT_CONNECT": ["api", "timeout_connect"],
    "CHEMBL_DA_TIMEOUT_READ": ["api", "timeout_read"],
    "CHEMBL_DA_OPENALEX_BASE": ["openalex", "base"],
    "CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT": ["openalex", "timeout_connect"],
    "CHEMBL_DA_OPENALEX_TIMEOUT_READ": ["openalex", "timeout_read"],
    "CHEMBL_DA_OPENALEX_RPS": ["openalex", "rps"],
    "CHEMBL_DA_OPENALEX_BURST": ["openalex", "burst"],
    "CHEMBL_DA_OPENALEX_MAILTO": ["openalex", "mailto"],
    "CHEMBL_DA_CROSSREF_BASE": ["crossref", "base"],
    "CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT": ["crossref", "timeout_connect"],
    "CHEMBL_DA_CROSSREF_TIMEOUT_READ": ["crossref", "timeout_read"],
    "CHEMBL_DA_CROSSREF_RPS": ["crossref", "rps"],
    "CHEMBL_DA_CROSSREF_BURST": ["crossref", "burst"],
    "CHEMBL_DA_CROSSREF_MAILTO": ["crossref", "mailto"],
    "CHEMBL_DA_UNIPROT_BASE": ["uniprot", "base"],
    "CHEMBL_DA_UNIPROT_TIMEOUT_CONNECT": ["uniprot", "timeout_connect"],
    "CHEMBL_DA_UNIPROT_TIMEOUT_READ": ["uniprot", "timeout_read"],
    "CHEMBL_DA_UNIPROT_RPS": ["uniprot", "rps"],
    "CHEMBL_DA_UNIPROT_BURST": ["uniprot", "burst"],
    "CHEMBL_DA_UNIPROT_MAPPING_BASE": ["uniprot_mapping", "base"],
    "CHEMBL_DA_UNIPROT_MAPPING_POLL_INTERVAL": ["uniprot_mapping", "poll_interval"],
    "CHEMBL_DA_UNIPROT_MAPPING_TIMEOUT": ["uniprot_mapping", "timeout"],
    "CHEMBL_DA_IUPHAR_BASE": ["iuphar", "base"],
    "CHEMBL_DA_IUPHAR_TIMEOUT_CONNECT": ["iuphar", "timeout_connect"],
    "CHEMBL_DA_IUPHAR_TIMEOUT_READ": ["iuphar", "timeout_read"],
    "CHEMBL_DA_IUPHAR_RPS": ["iuphar", "rps"],
    "CHEMBL_DA_IUPHAR_BURST": ["iuphar", "burst"],
    "CHEMBL_DA_PUBCHEM_BASE": ["pubchem", "base"],
    "CHEMBL_DA_PUBCHEM_TIMEOUT_CONNECT": ["pubchem", "timeout_connect"],
    "CHEMBL_DA_PUBCHEM_TIMEOUT_READ": ["pubchem", "timeout_read"],
    "CHEMBL_DA_PUBCHEM_RPS": ["pubchem", "rps"],
    "CHEMBL_DA_PUBCHEM_BURST": ["pubchem", "burst"],
    "CHEMBL_DA_OUTDIR": ["io", "output_dir"],
    "CHEMBL_DA_CACHE_DIR": ["io", "cache_dir"],
    "CHEMBL_DA_DICT_DIR": ["resources", "dictionary_dir"],
    "CHEMBL_DA_IUPHAR_TARGET_CSV": ["resources", "iuphar_target_csv"],
    "CHEMBL_DA_IUPHAR_FAMILY_CSV": ["resources", "iuphar_family_csv"],
    "CHEMBL_DA_UNIPROT_DATA_DIR": ["resources", "uniprot_data_dir"],
    "CHEMBL_DA_ORGANISM_CSV": ["resources", "organism_csv"],
    "CHEMBL_DA_STATUS_CSV": ["resources", "status_csv"],
    "CHEMBL_DA_TARGETS_TYPE_CSV": ["resources", "targets_type_csv"],
    "CHEMBL_DA_CONCURRENCY": ["jobs", "concurrency"],
    "CHEMBL_DA_CHUNK_SIZE": ["jobs", "chunk_size"],
    "CHEMBL_DA_GLOBAL_RPS": ["rate", "global_rps"],
    "CHEMBL_DA_GLOBAL_BURST": ["rate", "global_burst"],
    "CHEMBL_DA_LOG_LEVEL": ["log", "level"],
    "CHEMBL_DA_LOG_FORMAT": ["log", "format"],
    "CHEMBL_DA_LOG_DATEFMT": ["log", "datefmt"],
    "CHEMBL_DA_RETRY_MAX_ATTEMPTS": ["retry", "max_attempts"],
    "CHEMBL_DA_RETRY_BACKOFF_FACTOR": ["retry", "backoff_factor"],
}


def _apply_env_overrides(cfg: Config) -> None:
    """Apply environment variable overrides to ``cfg``.

    Environment variables that do not map to known configuration paths are
    ignored and generate a warning.
    """

    prefix = "CHEMBL_DA"
    for env_key, env_val in os.environ.items():
        key = env_key.upper()
        if key in _ALIAS_MAP:
            path = _ALIAS_MAP[key]
            try:
                _set_by_path(cfg, path, env_val)
            except KeyError:
                logger.warning(
                    "Environment variable %s ignored: unknown config path %s",
                    key,
                    ".".join(path),
                )
            continue
        if not key.startswith(prefix + "__"):
            continue
        parts = key.split("__")[1:]
        path_parts = [p.lower() for p in parts]
        try:
            _set_by_path(cfg, path_parts, env_val)
        except KeyError:
            logger.warning(
                "Environment variable %s ignored: unknown config path %s",
                key,
                ".".join(path_parts),
            )


def _serialize_paths(data: Any) -> Any:
    """Recursively convert :class:`~pathlib.Path` objects to strings."""

    if isinstance(data, dict):
        return {k: _serialize_paths(v) for k, v in data.items()}
    if isinstance(data, list):  # handle lists such as retry.status_forcelist
        return [_serialize_paths(v) for v in data]
    if isinstance(data, Path):
        return str(data)
    return data


def _mask_secrets(data: Any) -> Any:
    """Recursively mask values for keys that look like secrets.

    Keys containing common secret keywords (``key``, ``token``, ``secret``,
    ``password``) have their values replaced with ``"***"``. The check is
    case-insensitive and operates on nested mappings.
    """

    secret_tokens = {"key", "token", "secret", "password"}
    if isinstance(data, dict):
        masked: Dict[str, Any] = {}
        for k, v in data.items():
            if any(tok in k.lower() for tok in secret_tokens):
                masked[k] = "***"
            else:
                masked[k] = _mask_secrets(v)
        return masked
    if isinstance(data, list):
        return [_mask_secrets(v) for v in data]
    return data


def print_config(cfg: Config) -> None:
    """Print ``cfg`` in YAML format masking secret values.

    Parameters
    ----------
    cfg:
        Configuration object to serialise.
    """

    data = _serialize_paths(cfg.to_dict())
    masked = _mask_secrets(data)
    print(yaml.safe_dump(masked, sort_keys=False))


# JSON schema describing the configuration structure.
CONFIG_SCHEMA: Dict[str, Any] = {
    "type": "object",
    "properties": {
        "api": {
            "type": "object",
            "properties": {
                "chembl_base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "backoff_factor": {"type": "number", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
                "user_agent": {"type": "string"},
            },
            "required": [
                "chembl_base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "backoff_factor",
                "rps",
                "burst",
                "user_agent",
            ],
            "additionalProperties": False,
        },
        "openalex": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
                "mailto": {"type": "string"},
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "rps",
                "burst",
                "mailto",
            ],
            "additionalProperties": False,
        },
        "crossref": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
                "mailto": {"type": "string"},
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "rps",
                "burst",
                "mailto",
            ],
            "additionalProperties": False,
        },
        "uniprot": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
                "delay": {"type": "number", "minimum": 0},
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "rps",
                "burst",
            ],
            "additionalProperties": False,
        },
        "uniprot_mapping": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "poll_interval": {"type": "number", "exclusiveMinimum": 0},
                "timeout": {"type": "number", "minimum": 1},
            },
            "required": ["base", "poll_interval", "timeout"],
            "additionalProperties": False,
        },
        "iuphar": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "rps",
                "burst",
            ],
            "additionalProperties": False,
        },
        "pubchem": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "rps": {"type": "integer", "minimum": 1},
                "burst": {"type": "integer", "minimum": 1},
                "delay": {"type": "number", "minimum": 0},
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "rps",
                "burst",
            ],
            "additionalProperties": False,
        },
        "pubmed": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "encodings": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                },
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "encodings",
            ],
            "additionalProperties": False,
        },
        "semantic_scholar": {
            "type": "object",
            "properties": {
                "base": {"type": "string", "format": "uri"},
                "timeout_connect": {"type": "integer", "minimum": 1},
                "timeout_read": {"type": "integer", "minimum": 1},
                "retries": {"type": "integer", "minimum": 0},
                "encodings": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                },
            },
            "required": [
                "base",
                "timeout_connect",
                "timeout_read",
                "retries",
                "encodings",
            ],
            "additionalProperties": False,
        },
        "doc_type": {
            "type": "object",
            "properties": {
                "weights": {
                    "type": "object",
                    "properties": {
                        "pubmed": {"type": "integer", "minimum": 0},
                        "openalex": {"type": "integer", "minimum": 0},
                        "scholar": {"type": "integer", "minimum": 0},
                    },
                    "required": ["pubmed", "openalex", "scholar"],
                    "additionalProperties": False,
                },
                "thresholds": {
                    "type": "object",
                    "properties": {
                        "review": {"type": "integer", "minimum": 0},
                        "experimental": {"type": "integer", "minimum": 0},
                        "unknown": {"type": "integer", "minimum": 0},
                    },
                    "required": ["review", "experimental", "unknown"],
                    "additionalProperties": False,
                },
            },
            "required": ["weights", "thresholds"],
            "additionalProperties": False,
        },
        "resources": {
            "type": "object",
            "properties": {
                "dictionary_dir": {"type": "string", "minLength": 1},
                "iuphar_target_csv": {"type": "string", "minLength": 1},
                "iuphar_family_csv": {"type": "string", "minLength": 1},
                "uniprot_data_dir": {"type": "string", "minLength": 1},
                "organism_csv": {"type": "string", "minLength": 1},
                "status_csv": {"type": "string", "minLength": 1},
                "targets_type_csv": {"type": "string", "minLength": 1},
            },
            "required": [
                "dictionary_dir",
                "iuphar_target_csv",
                "iuphar_family_csv",
                "uniprot_data_dir",
                "organism_csv",
                "status_csv",
                "targets_type_csv",
            ],
            "additionalProperties": False,
        },
        "io": {
            "type": "object",
            "properties": {
                "output_dir": {"type": "string", "minLength": 1},
                "cache_dir": {"type": "string", "minLength": 1},
                "csv_sep": {"type": "string", "minLength": 1},
                "csv_encoding": {"type": "string", "minLength": 1},
                "csv_index": {"type": "boolean"},
                "exist_ok": {"type": "boolean"},
            },
            "required": [
                "output_dir",
                "cache_dir",
                "csv_sep",
                "csv_encoding",
                "csv_index",
                "exist_ok",
            ],
            "additionalProperties": False,
        },
        "jobs": {
            "type": "object",
            "properties": {
                "concurrency": {"type": "integer", "minimum": 1},
                "chunk_size": {"type": "integer", "minimum": 1},
            },
            "required": ["concurrency", "chunk_size"],
            "additionalProperties": False,
        },
        "batch": {
            "type": "object",
            "properties": {
                "size": {"type": "integer", "minimum": 1},
                "pause": {"type": "number", "minimum": 0},
                "concurrency": {"type": "integer", "minimum": 1},
                "fail_fast": {"type": "boolean"},
                "retry_failed": {"type": "boolean"},
            },
            "required": [
                "size",
                "pause",
                "concurrency",
                "fail_fast",
                "retry_failed",
            ],
            "additionalProperties": False,
        },
        "quality": {
            "type": "object",
            "properties": {
                "sample_rows": {"type": "integer", "minimum": 0},
                "corr_method": {"type": "string", "minLength": 1},
                "max_unique_preview": {"type": "integer", "minimum": 1},
                "bin_count": {"type": "integer", "minimum": 1},
            },
            "required": [
                "sample_rows",
                "corr_method",
                "max_unique_preview",
                "bin_count",
            ],
            "additionalProperties": False,
        },
        "mapper": {
            "type": "object",
            "properties": {
                "enable_cache": {"type": "boolean"},
                "strict_schema": {"type": "boolean"},
                "warn_on_cast": {"type": "boolean"},
            },
            "required": ["enable_cache", "strict_schema", "warn_on_cast"],
            "additionalProperties": False,
        },
        "init": {
            "type": "object",
            "properties": {
                "same_doc": {"type": "string", "minLength": 1},
                "all_doc": {"type": "string", "minLength": 1},
                "output_dir": {"type": "string", "minLength": 1},
            },
            "required": ["same_doc", "all_doc", "output_dir"],
            "additionalProperties": False,
        },
        "rate": {
            "type": "object",
            "properties": {
                "global_rps": {"type": "integer", "minimum": 1},
                "global_burst": {"type": "integer", "minimum": 1},
            },
            "required": ["global_rps", "global_burst"],
            "additionalProperties": False,
        },
        "retry": {
            "type": "object",
            "properties": {
                "max_attempts": {"type": "integer", "minimum": 1},
                "backoff_factor": {"type": "number", "minimum": 0},
                "status_forcelist": {
                    "type": "array",
                    "items": {"type": "integer"},
                },
            },
            "required": ["max_attempts", "backoff_factor", "status_forcelist"],
            "additionalProperties": False,
        },
        "log": {
            "type": "object",
            "properties": {
                "level": {"type": "string", "minLength": 1},
                "format": {"type": "string", "minLength": 1},
                "datefmt": {"type": "string", "minLength": 1},
            },
            "required": ["level", "format", "datefmt"],
            "additionalProperties": False,
        },
    },
    "required": [
        "api",
        "openalex",
        "crossref",
        "uniprot",
        "uniprot_mapping",
        "iuphar",
        "pubchem",
        "pubmed",
        "semantic_scholar",
        "doc_type",
        "resources",
        "io",
        "jobs",
        "batch",
        "quality",
        "mapper",
        "init",
        "rate",
        "retry",
        "log",
    ],
    "additionalProperties": False,
}


def _validate(cfg: Config) -> None:
    """Validate ``cfg`` against :data:`CONFIG_SCHEMA`."""

    validator = jsonschema.Draft202012Validator(
        CONFIG_SCHEMA, format_checker=jsonschema.FormatChecker()
    )
    # ``cfg`` contains ``Path`` instances; convert them to strings before validation
    validator.validate(_serialize_paths(cfg.to_dict()))

    # Validate logging level (case-insensitive)
    # ``logging.getLevelNamesMapping`` was added in Python 3.11. Fallback to the
    # private ``logging._nameToLevel`` mapping for older versions.
    try:
        level_names = logging.getLevelNamesMapping()
    except AttributeError:  # pragma: no cover - python <3.11 only
        level_names = {
            name.upper(): level for name, level in logging._nameToLevel.items()
        }
    if cfg.log.level.upper() not in level_names:
        valid = ", ".join(sorted(level_names))
        raise ValueError(f"log.level must be one of {valid}, got {cfg.log.level!r}")

    """Basic sanity checks for configuration values."""
    if not _valid_url(cfg.api.chembl_base):
        raise ValueError("api.chembl_base must be a valid URL")
    if cfg.api.timeout_connect <= 0 or cfg.api.timeout_read <= 0:
        raise ValueError("api timeouts must be positive")
    if cfg.api.retries < 0 or cfg.api.backoff_factor < 0:
        raise ValueError(
            "api.retries must be non-negative and backoff_factor non-negative"
        )
    if cfg.api.rps <= 0 or cfg.api.burst <= 0:
        raise ValueError("api.rps and api.burst must be positive")
    if not _EMAIL_RE.search(cfg.api.user_agent):
        raise ValueError(
            "api.user_agent must include contact information such as an email"
        )

    services_full: list[tuple[str, Any]] = [
        ("openalex", cfg.openalex),
        ("crossref", cfg.crossref),
        ("uniprot", cfg.uniprot),
        ("iuphar", cfg.iuphar),
        ("pubchem", cfg.pubchem),
    ]
    for name, service in services_full:
        if not _valid_url(service.base):
            raise ValueError(f"{name}.base must be a valid URL")
        if service.timeout_connect <= 0 or service.timeout_read <= 0:
            raise ValueError(f"{name} timeouts must be positive")
        if service.retries < 0:
            raise ValueError(f"{name}.retries must be non-negative")
        if service.rps <= 0 or service.burst <= 0:
            raise ValueError(f"{name}.rps and {name}.burst must be positive")
        if hasattr(service, "delay") and service.delay < 0:
            raise ValueError(f"{name}.delay must be non-negative")

    basic_services: list[tuple[str, Any]] = [
        ("pubmed", cfg.pubmed),
        ("semantic_scholar", cfg.semantic_scholar),
    ]
    for name, service in basic_services:
        if not _valid_url(service.base):
            raise ValueError(f"{name}.base must be a valid URL")
        if service.timeout_connect <= 0 or service.timeout_read <= 0:
            raise ValueError(f"{name} timeouts must be positive")
        if service.retries < 0:
            raise ValueError(f"{name}.retries must be non-negative")
        if not service.encodings:
            raise ValueError(f"{name}.encodings must not be empty")

    mapping = cfg.uniprot_mapping
    if not _valid_url(mapping.base):
        raise ValueError("uniprot_mapping.base must be a valid URL")
    if mapping.poll_interval <= 0 or mapping.timeout <= 0:
        raise ValueError("uniprot_mapping.poll_interval and timeout must be positive")

    for name, mail in [
        ("openalex", cfg.openalex.mailto),
        ("crossref", cfg.crossref.mailto),
    ]:
        if not mail or not _EMAIL_RE.fullmatch(mail):
            raise ValueError(f"{name}.mailto must be a valid email address")

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

    out_dir = cfg.io.output_dir
    cache_dir = cfg.io.cache_dir
    for path in (out_dir, cache_dir):
        if not path.exists() and not cfg.io.exist_ok:
            raise FileNotFoundError(f"{path} does not exist")
        if path.exists() and not path.is_dir():
            raise NotADirectoryError(f"{path} is not a directory")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def ensure_dirs(cfg: Config) -> None:
    """Create I/O directories if required.

    Parameters
    ----------
    cfg : Config
        Configuration settings containing ``io`` paths.

    Raises
    ------
    FileNotFoundError
        If a directory is missing and ``cfg.io.exist_ok`` is ``False``.
    NotADirectoryError
        If an existing path is not a directory.
    """

    out_dir = cfg.io.output_dir
    cache_dir = cfg.io.cache_dir
    for path in (out_dir, cache_dir):
        if path.exists():
            if not path.is_dir():
                raise NotADirectoryError(f"{path} is not a directory")
        else:
            if cfg.io.exist_ok:
                path.mkdir(parents=True, exist_ok=True)
            else:
                raise FileNotFoundError(f"{path} does not exist")


def load_config(
    path: str | Path = "config.yaml",
    cli_overrides: Dict[str, Any] | None = None,
    *,
    strict: bool = False,
) -> Config:
    """Load configuration from ``path`` applying environment and CLI overrides.

    Parameters
    ----------
    path:
        Location of the YAML configuration file.
    cli_overrides:
        Mapping of ``"section.key"`` paths to values coming from the command
        line.
    strict:
        When ``True`` raise :class:`ValueError` for unknown configuration keys;
        otherwise a warning is emitted.

    Returns
    -------
    Config
        Fully populated configuration object.

    Raises
    ------
    ConfigError
        If ``path`` does not exist or contains invalid YAML.
    ValueError
        If ``strict`` is ``True`` and unknown configuration keys are present.
    TypeError
        If configuration values have incorrect types.
    """

    cfg = Config()
    unknown_keys: List[str] = []
    try:
        with Path(path).open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}
    except FileNotFoundError as err:
        raise ConfigError(f"configuration file not found: {path}") from err
    except yaml.YAMLError as err:
        raise ConfigError(
            f"failed to parse YAML configuration at {path}: {err}"
        ) from err
    if not isinstance(data, dict):
        raise TypeError("top-level structure in config file must be a mapping")
    _update_from_dict(cfg, data, unknown_keys=unknown_keys)

    if unknown_keys:
        joined = ", ".join(sorted(unknown_keys))
        msg = f"Unknown configuration key(s) in {path}: {joined}"
        if strict:
            raise ValueError(msg)
        logger.warning(msg)

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
    "UniprotMappingCfg",
    "IupharCfg",
    "PubChemCfg",
    "PubMedCfg",
    "SemanticScholarCfg",
    "DocTypeCfg",
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
]
