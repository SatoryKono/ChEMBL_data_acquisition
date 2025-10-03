"""Application configuration using Pydantic models.

This module provides a typed wrapper around ``config/config.yaml``. Configuration
values are loaded from a YAML file and can be overridden by environment
variables or command line options. The order of precedence is::

    YAML < config.local.yaml < environment variables < CLI overrides

Environment variables follow the ``CHEMBL_DA__SECTION__KEY`` naming pattern
where sections and keys are joined by double underscores. A number of short
aliases are supported; see ``_ALIAS_MAP`` for the full list.
"""

from __future__ import annotations

import atexit
import logging
import os
import re
from collections.abc import Sequence
from contextlib import ExitStack
from importlib import resources
from pathlib import Path
from types import UnionType
from typing import Any, Mapping, Union, get_args, get_origin
from urllib.parse import urlparse

import yaml
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    StrictInt,
    ValidationError,
    field_validator,
)
from pydantic_core import ErrorDetails
from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .common.rate_limiter import configure_limiter_cache

from .common.log import logger
from .utils.config import ConfigLoaderError, load_yaml_config


def _annotation_is_path(annotation: Any) -> bool:
    """Return ``True`` if *annotation* represents a :class:`~pathlib.Path`."""

    try:
        return issubclass(annotation, Path)
    except TypeError:
        return False


def _annotation_is_model(annotation: Any) -> bool:
    """Return ``True`` if *annotation* is a :class:`pydantic.BaseModel`."""

    try:
        return issubclass(annotation, BaseModel)
    except TypeError:
        return False


def _collect_path_field_paths(
    model: type[BaseModel], prefix: tuple[str, ...] = ()
) -> set[tuple[str, ...]]:
    """Return dotted paths for ``Path`` fields within *model*."""

    paths: set[tuple[str, ...]] = set()
    for name, field in model.model_fields.items():
        annotation = field.annotation
        if annotation is None:
            continue
        origin = get_origin(annotation)
        if origin in {UnionType, Union}:
            args = [arg for arg in get_args(annotation) if arg is not type(None)]
            if any(_annotation_is_path(arg) for arg in args):
                paths.add((*prefix, name))
                continue
            if any(_annotation_is_model(arg) for arg in args):
                for arg in args:
                    if _annotation_is_model(arg):
                        paths.update(_collect_path_field_paths(arg, (*prefix, name)))
                continue
        if _annotation_is_path(annotation):
            paths.add((*prefix, name))
            continue
        if _annotation_is_model(annotation):
            paths.update(_collect_path_field_paths(annotation, (*prefix, name)))
    return paths


def _absolutise_path_value(value: Any, base_dir: Path) -> Any:
    """Return *value* converted to an absolute path relative to *base_dir*."""

    if value is None:
        return value
    if isinstance(value, str):
        path = Path(value)
        if path.is_absolute():
            return value
        return str((base_dir / path).resolve())
    if isinstance(value, os.PathLike):
        path = Path(value)
        if path.is_absolute():
            return path
        return (base_dir / path).resolve()
    return value


def _absolutise_config_paths(data: Mapping[str, Any], base_dir: Path) -> None:
    """Normalise relative ``Path`` entries in *data* using *base_dir*."""

    for path in _CONFIG_PATH_FIELDS:
        current: Mapping[str, Any] | Any = data
        for key in path[:-1]:
            if not isinstance(current, Mapping):
                current = None
                break
            current = current.get(key)
        if not isinstance(current, Mapping):
            continue
        final_key = path[-1]
        if final_key not in current:
            continue
        value = current[final_key]
        if isinstance(current, dict):
            current[final_key] = _absolutise_path_value(value, base_dir)


def _normalize_base_path(value: Path | str) -> Path:
    """Return *value* coerced into an absolute :class:`~pathlib.Path`."""

    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    return (Path.cwd() / candidate).resolve()


def _default_base_path() -> Path:
    """Return the default data directory used for placeholder expansion."""

    env_override = os.environ.get("CHEMBL_DA_BASE_PATH")
    if env_override:
        return _normalize_base_path(env_override)
    return (Path.home() / ".local" / "share" / "chembl-da").resolve()


def _default_cache_home() -> Path:
    """Return the default cache directory for local artefacts."""

    return (Path.home() / ".cache" / "chembl-da").resolve()


def _expand_config_placeholders(data: Any, *, base_path: Path) -> Any:
    """Expand configuration placeholders in ``data`` using *base_path*."""

    replacements = {"$CHEMBL_DA_BASE_PATH": str(base_path)}

    def _expand(value: Any) -> Any:
        if isinstance(value, dict):
            for key, current in value.items():
                value[key] = _expand(current)
            return value
        if isinstance(value, list):
            for idx, current in enumerate(value):
                value[idx] = _expand(current)
            return value
        if isinstance(value, tuple):
            return tuple(_expand(item) for item in value)
        if isinstance(value, str):
            replaced = value
            for marker, target in replacements.items():
                replaced = replaced.replace(marker, target)
            replaced = os.path.expandvars(replaced)
            replaced = os.path.expanduser(replaced)
            return replaced
        return value

    return _expand(data)


def _resolve_placeholder_base_path(base_path: Path | str | None) -> Path:
    """Return the base path used for ``$CHEMBL_DA_BASE_PATH`` placeholders."""

    if base_path is not None:
        return _normalize_base_path(base_path)
    return _default_base_path()


_RESOURCE_STACK = ExitStack()
atexit.register(_RESOURCE_STACK.close)


def _dictionary_resource(*parts: str) -> Path:
    """Return a filesystem path for a bundled dictionary resource."""

    traversable = resources.files("dictionary")
    for part in parts:
        traversable = traversable.joinpath(part)
    return Path(_RESOURCE_STACK.enter_context(resources.as_file(traversable)))


_EMAIL_RE = re.compile(r"[^@\s]+@[^@\s]+\.[^@\s]+")
_PLACEHOLDER_EMAIL_DOMAINS = {
    "example.com",
    "example.net",
    "example.org",
    "test.com",
    "test.net",
    "test.org",
}


def _require_non_placeholder_email(service: str, value: str) -> str:
    """Ensure *value* is a real contact email for *service*."""

    if not value or not _EMAIL_RE.fullmatch(value):
        raise ValueError(f"{service}.mailto must be a valid email")
    _, domain = value.rsplit("@", 1)
    domain = domain.lower()
    if domain in _PLACEHOLDER_EMAIL_DOMAINS or domain.endswith(".example"):
        raise ValueError(
            f"{service}.mailto must not use placeholder domain {domain}; "
            "configure a real contact email."
        )
    return value


def _valid_url(url: str) -> bool:
    """Return ``True`` if *url* contains both scheme and netloc components."""

    parsed = urlparse(url)
    return bool(parsed.scheme and parsed.netloc)


def _coerce_integral_numbers(value: Any) -> Any:
    """Return *value* with floats that represent integers converted to ``int``."""

    if isinstance(value, float):
        return int(value) if value.is_integer() else value
    if isinstance(value, list):
        return [_coerce_integral_numbers(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_coerce_integral_numbers(item) for item in value)
    if isinstance(value, dict):
        return {key: _coerce_integral_numbers(val) for key, val in value.items()}
    return value


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


# ChEMBL and PubChem test item fields requested by default when fetching molecule data.
TESTITEM_FIELD_DEFAULTS: tuple[str, ...] = (
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "pref_name",
    "max_phase",
    "molecule_type",
    "first_approval",
    "oral",
    "parenteral",
    "topical",
    "black_box_warning",
    "structure_type",
    "molecule_structures.canonical_smiles",
    "molecule_structures.standard_inchi",
    "molecule_structures.standard_inchi_key",
    "pubchem_cid",
    "pubchem_iupac_name",
    "pubchem_molecular_formula",
    "pubchem_isomeric_smiles",
    "pubchem_canonical_smiles",
    "pubchem_inchi",
    "pubchem_inchikey",
)


class ApiCfg(_BaseModel):
    """Settings for ChEMBL API access."""

    chembl_base: str = "https://www.ebi.ac.uk/chembl/api/data"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(30.0, ge=1)
    retries: int = Field(3, ge=0)
    backoff_factor: float = Field(0.5, ge=0)
    rps: int = Field(5, ge=1)
    burst: int = Field(5, ge=1)
    user_agent: str = "chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)"

    @field_validator("chembl_base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):  # pragma: no cover - defensive
            raise ValueError("invalid URL")
        return v

    @field_validator("user_agent")
    @classmethod
    def _ua(cls, v: str) -> str:
        if "contact@example.org" in v:
            raise ConfigError(
                "api.user_agent must include a real contact address; "
                "replace contact@example.org with a valid email"
            )
        if not _EMAIL_RE.search(v):
            raise ValueError(
                "api.user_agent must include contact information such as an email"
            )
        return v


class ChemblCacheCfg(_BaseModel):
    cache_ttl: int = Field(3600, ge=1)
    cache_maxsize: int = Field(1024, ge=1)


class MoleculeCatalogCfg(_BaseModel):
    cache_path: Path = Field(
        default_factory=lambda: _default_base_path()
        / "cache"
        / "molecule_parent_catalog.json"
    )
    sqlite_path: Path = Field(
        default_factory=lambda: _default_base_path()
        / "cache"
        / "molecule_parent_catalog.sqlite"
    )
    endpoint: str = "molecule"
    child_field: str = "molecule_chembl_id"
    parent_field: str = "parent_molecule_chembl_id"
    hierarchy_lookup_path: Path | None = Path(
        "dictionary/_testitem/molecule_hierarchy.csv"
    )
    hierarchy_lookup_encoding: str = "utf-8-sig"
    hierarchy_lookup_delimiter: str = ","
    force_refresh_existing: bool = False
    fields: tuple[str, ...] = (
        "molecule_chembl_id",
        "parent_molecule_chembl_id",
    )
    filters: dict[str, str] = Field(
        default_factory=lambda: {"parent_molecule_chembl_id__isnull": "false"}
    )
    page_size: int = Field(500, ge=1)
    fallback_single_limit: int | None = Field(default=None, ge=1)


class OpenAlexCfg(_BaseModel):
    base: str = "https://api.openalex.org"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(30.0, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)
    mailto: str = "chembl-data@ebi.ac.uk"

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    @field_validator("mailto")
    @classmethod
    def _email(cls, v: str) -> str:
        return _require_non_placeholder_email("openalex", v)


class CrossRefCfg(_BaseModel):
    base: str = "https://api.crossref.org"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(30.0, ge=1)
    retries: int = Field(3, ge=0)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)
    mailto: str = "chembl-data@ebi.ac.uk"

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    @field_validator("mailto")
    @classmethod
    def _email(cls, v: str) -> str:
        return _require_non_placeholder_email("crossref", v)


class UniprotCfg(_BaseModel):
    base: str = "https://rest.uniprot.org"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(30.0, ge=1)
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
                self.cache_ttl,
            )
        )


class IupharCfg(_BaseModel):
    base: str = "https://www.guidetopharmacology.org/services"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(30.0, ge=1)
    rps: int = Field(4, ge=1)
    burst: int = Field(5, ge=1)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class PubChemCfg(_BaseModel):
    """Settings for resolving PubChem data."""

    enable: bool = Field(
        True,
        description="Enable PubChem augmentation for test item data",
    )
    base: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    user_agent: str = Field(
        "chembl-da/1.0 (mailto:chembl-data@ebi.ac.uk)",
        description="Custom User-Agent for PubChem requests including contact details",
    )
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(60.0, ge=1)
    timeout_seconds: float = Field(
        30.0,
        ge=0,
        description="Overall timeout applied to individual PubChem lookups",
    )
    retries: int = Field(3, ge=0)
    rps: int = Field(3, ge=1)
    burst: int = Field(5, ge=1)
    delay: float = Field(
        0.2,
        ge=0,
        description="Base delay between retries for network errors",
    )
    backoff_initial_seconds: float = Field(
        0.5,
        ge=0,
        description="Initial delay when backing off after 429/5xx responses",
    )
    resolve_order: tuple[str, ...] = Field(
        (
            "cache",
            "smiles",
            "inchikey",
            "inchi",
            "pref_name",
        ),
        description="Ordered list of lookup strategies when resolving PubChem CIDs",
    )
    cache_ttl: int = Field(
        3600,
        ge=0,
        description="Time-to-live for the in-memory PubChem request cache in seconds",
    )
    cache_maxsize: int = Field(
        1024,
        ge=1,
        description="Maximum number of entries stored in the PubChem in-memory cache",
    )
    cache_ttl_hours: float | None = Field(
        None,
        ge=0,
        description="Optional TTL for the persisted CID cache in hours",
    )
    cid_cache_path: Path | None = Field(
        default_factory=lambda: _default_base_path()
        / "cache"
        / "pubchem_cid_cache.json",
        description="Optional JSON cache storing PubChem CIDs by molecule_chembl_id",
    )
    batch_size: int = Field(
        50,
        ge=1,
        description="Maximum PubChem property lookups submitted per batch; concurrency is capped by pubchem.rps",
    )
    prefer_local_smiles: bool = Field(
        False,
        description="Skip PubChem lookups when existing pubchem_* columns are complete",
    )
    prefer_local_values: bool = Field(
        True,
        description="Retain existing pubchem_* values when lookups return empty results",
    )
    use_parent_for_salts: bool = Field(
        True,
        description="Resolve PubChem CIDs via parent structures when child lookups fail",
    )
    allow_polymer: bool = Field(
        False,
        description="Allow PubChem lookups for polymer or mixture records",
    )
    write_not_found_literal: bool = Field(
        False,
        description="Write 'Not Found' when PubChem lookups fail",
    )

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v

    @field_validator("user_agent")
    @classmethod
    def _ua(cls, v: str) -> str:
        if not _EMAIL_RE.search(v):
            raise ValueError(
                "pubchem.user_agent must include contact information such as an email",
            )
        return v


class PubMedCfg(_BaseModel):
    base: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(10.0, ge=1)
    retries: int = Field(2, ge=0)
    rps: int | None = Field(default=None, ge=1)
    burst: int | None = Field(default=None, ge=1)

    @field_validator("base")
    @classmethod
    def _url(cls, v: str) -> str:
        if not _valid_url(v):
            raise ValueError("invalid URL")
        return v


class SemanticScholarCfg(_BaseModel):
    base: str = "https://api.semanticscholar.org/graph/v1"
    timeout_connect: float = Field(5.0, ge=1)
    timeout_read: float = Field(10.0, ge=1)
    retries: int = Field(2, ge=0)
    rps: int | None = Field(default=None, ge=1)
    burst: int | None = Field(default=None, ge=1)

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
    limit: int | None = Field(default=None, ge=0)


class DocQualityCfg(_BaseModel):
    enable: bool = True
    sample_rows: int | None = Field(default=None, ge=1)
    include_columns: tuple[str, ...] | None = None
    exclude_columns: tuple[str, ...] | None = None
    fatal_on_error: bool = False


class ResourcesCfg(_BaseModel):
    dictionary_dir: Path = Field(default_factory=_dictionary_resource)
    iuphar_target_csv: Path = Field(
        default_factory=lambda: _dictionary_resource(
            "_target", "_IUPHAR", "_IUPHAR_target.csv"
        )
    )
    iuphar_family_csv: Path = Field(
        default_factory=lambda: _dictionary_resource(
            "_target", "_IUPHAR", "_IUPHAR_family.csv"
        )
    )
    uniprot_data_dir: Path = Field(
        default_factory=lambda: _dictionary_resource("_target", "_uniprot")
    )
    targets_type_csv: Path = Field(
        default_factory=lambda: _dictionary_resource("_target", "targets_type.csv")
    )

    @field_validator(
        "dictionary_dir",
        "iuphar_target_csv",
        "iuphar_family_csv",
        "uniprot_data_dir",
        "targets_type_csv",
        mode="before",
    )
    @classmethod
    def _coerce_path(cls, value: Any) -> Path | Any:
        if isinstance(value, (str, os.PathLike)):
            return Path(value)
        return value


class IoCfg(_BoolModel):
    output_dir: Path = Field(
        default_factory=lambda: _default_base_path() / "output"
    )
    cache_dir: Path = Field(
        default_factory=lambda: _default_cache_home()
    )
    csv_sep: str = ","
    csv_encoding: str = "utf-8-sig"
    csv_fallback_encodings: Sequence[str] | None = Field(
        default_factory=lambda: ("utf-8", "cp1252", "windows-1251", "latin-1")
    )
    csv_fallback_separators: Sequence[str] | None = Field(
        default_factory=lambda: ("\t", ";")
    )
    na_markers: Sequence[str] | None = ("#N/A",)
    csv_chunksize: int = Field(10000, ge=1)
    exist_ok: bool = True
    keep_na_markers: bool = False

    @field_validator("exist_ok", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class LogCfg(_BaseModel):
    level: str = "INFO"

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
    same_doc: Path = Field(
        default_factory=lambda: _default_base_path()
        / "input"
        / "ChEMBL"
        / "ChEMBL_same_document_20_05.xlsx"
    )
    all_doc: Path = Field(
        default_factory=lambda: _default_base_path()
        / "input"
        / "ChEMBL"
        / "ChEMBL_all_10_05_step5.xlsx"
    )
    output_dir: Path = Field(
        default_factory=lambda: _default_base_path()
        / "output"
        / "ChEMBL"
        / "processed"
    )


class RateCfg(_BaseModel):
    global_rps: int | None = Field(8, ge=0)
    global_burst: int | None = Field(8, ge=0)
    limiter_cache_maxsize: int = Field(128, ge=1)
    limiter_cache_ttl: int = Field(600, ge=1)

    @field_validator("global_rps", "global_burst", mode="before")
    @classmethod
    def _allow_zero_or_none(cls, value: Any) -> Any:
        """Normalise disabled limiter values."""

        if value is None:
            return None
        if isinstance(value, str):
            stripped = value.strip()
            if not stripped:
                return None
            lowered = stripped.lower()
            if lowered in {"none", "null"}:
                return None
        return value


class RetryCfg(_BaseModel):
    max_attempts: int = Field(3, ge=1)
    backoff_factor: float = Field(0.5, ge=0)
    backoff_cap: float | None = Field(default=None, ge=0)
    status_forcelist: list[StrictInt] = Field(
        default_factory=lambda: [429, 500, 502, 503, 504]
    )


class ActivityBoundsCfg(_BoolModel):
    enable_from_relation: bool = True
    enable_from_uncertainty: bool = False
    rounding_digits: int = Field(3, ge=0)
    clamp_nonnegative: bool = True
    log_unknown_relations: bool = True

    @field_validator(
        "enable_from_relation",
        "enable_from_uncertainty",
        "clamp_nonnegative",
        "log_unknown_relations",
        mode="before",
    )
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class ActivityActionTypeCfg(_BoolModel):
    enabled: bool = True
    column: str = "action_type"
    log_missing: bool = True
    log_distribution: bool = True
    metrics: dict[str, str] = Field(
        default_factory=lambda: {
            "ic50": "inhibition",
            "ec50": "activation",
            "ac50": "activation",
            "ki": "binding",
            "kd": "binding",
        }
    )
    triages: dict[str, str] = Field(default_factory=dict)
    functionality: dict[str, str] = Field(
        default_factory=lambda: {
            "agonist": "activation",
            "partial agonist": "activation",
            "antagonist": "inhibition",
            "inhibitor": "inhibition",
            "activator": "activation",
        }
    )
    mechanism: dict[str, str] = Field(default_factory=dict)
    triage_fields: list[str] = Field(
        default_factory=lambda: [
            "data_validity_description",
            "data_validity_comment",
        ]
    )
    functionality_fields: list[str] = Field(
        default_factory=lambda: [
            "functional_activity",
            "action_type",
        ]
    )
    mechanism_fields: list[str] = Field(
        default_factory=lambda: [
            "mechanism_of_action",
            "mechanism_comment",
        ]
    )
    allowlist: list[str] = Field(
        default_factory=lambda: [
            "activation",
            "inhibition",
            "binding",
            "pam",
            "nam",
            "triaged",
            "unknown",
        ]
    )
    positive_label: str = "PAM"
    negative_label: str = "NAM"
    fallback: str | None = "unknown"

    @field_validator(
        "enabled",
        "log_missing",
        "log_distribution",
        mode="before",
    )
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)

    @field_validator("column")
    @classmethod
    def _column(cls, v: str) -> str:
        if not v or not v.strip():
            raise ValueError("activity_enrichment.action_type.column must be non-empty")
        return v

    @field_validator(
        "triage_fields",
        "functionality_fields",
        "mechanism_fields",
        mode="before",
    )
    @classmethod
    def _list(cls, value: Any) -> list[str]:
        if isinstance(value, str):
            items = [value]
        else:
            items = list(value)
        cleaned = [str(item).strip() for item in items if str(item).strip()]
        if not cleaned:
            raise ValueError(
                "activity_enrichment.action_type field lists must be non-empty"
            )
        return cleaned


class ActivityPropertiesCfg(_BoolModel):
    enabled: bool = True
    column: str = "activity_properties"
    summary_column: str = "activity_property_summary"
    name_field: str = "type"
    value_field: str = "value"
    units_field: str | None = "units"
    separator: str = "; "
    pair_separator: str = "="
    drop_source_column: bool = True
    log_missing: bool = False
    log_distribution: bool = False
    allowlist: list[str] = Field(
        default_factory=lambda: [
            "measurement",
            "assay",
            "comments",
            "effect_features",
            "triage",
            "mechanism",
            "functionality",
        ]
    )
    hash_column: str | None = "properties_hash"

    @field_validator(
        "enabled",
        "drop_source_column",
        "log_missing",
        "log_distribution",
        mode="before",
    )
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)

    @field_validator(
        "column",
        "summary_column",
        "name_field",
        "value_field",
        "separator",
        "pair_separator",
    )
    @classmethod
    def _non_empty(cls, v: str) -> str:
        if not v or not str(v).strip():
            raise ValueError(
                "activity_enrichment.activity_properties fields must be non-empty"
            )
        return v


class ActivityEnrichmentCfg(_BaseModel):
    action_type: ActivityActionTypeCfg = Field(
        default_factory=lambda: ActivityActionTypeCfg()
    )
    activity_properties: ActivityPropertiesCfg = Field(
        default_factory=lambda: ActivityPropertiesCfg()
    )


class ActivityCfg(_BoolModel):
    column: str = "activity_id"
    batch_size: int = Field(5, ge=1)
    workers: int = Field(1, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    dry_run: bool = False

    @field_validator("dry_run", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class AssayCfg(_BaseModel):
    column: str = "assay_chembl_id"
    batch_size: int = Field(10, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)


class TestitemBatchRetryCfg(_BoolModel):
    enable: bool = False
    shrink_factor: float = Field(0.5, gt=0, lt=1)
    min_size: int = Field(1, ge=1)

    @field_validator("enable", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class TestitemCfg(_BaseModel):
    column: str = "molecule_chembl_id"
    batch_size: int = Field(1000, ge=1, le=1000)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    fields: tuple[str, ...] = Field(default=TESTITEM_FIELD_DEFAULTS)
    request_limit: int = Field(1000, ge=1, le=1000)
    retries: int | None = Field(default=None, ge=0)
    backoff_factor: float | None = Field(default=None, ge=0)
    batch_retry: TestitemBatchRetryCfg = Field(
        default_factory=lambda: TestitemBatchRetryCfg()
    )

    @field_validator("fields", mode="before")
    @classmethod
    def _coerce_fields(cls, value: Any) -> tuple[str, ...]:
        if value is None:
            return TESTITEM_FIELD_DEFAULTS
        if isinstance(value, (list, tuple)):
            return tuple(str(item) for item in value)
        raise TypeError("testitem.fields must be a list or tuple of strings")


class TestitemMoleculeEnrichmentSourcesCfg(_BaseModel):
    molecule_catalog_path: Path = Path("dictionary/molecule_catalog.csv")
    molecule_hierarchy_path: Path = Path("dictionary/molecule_hierarchy.csv")


class TestitemMoleculeEnrichmentOutputCfg(_BoolModel):
    salt_as_null_when_absent: bool = True

    @field_validator("salt_as_null_when_absent", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class TestitemMoleculeEnrichmentFlagsCfg(_BoolModel):
    coerce_to_bool: bool = True
    parent_fallback: bool = True

    @field_validator("coerce_to_bool", "parent_fallback", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class TestitemMoleculeEnrichmentLoggingCfg(_BoolModel):
    warn_missing_molecule: bool = True
    warn_inconsistent_flags: bool = True

    @field_validator("warn_missing_molecule", "warn_inconsistent_flags", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class TestitemMoleculeEnrichmentCfg(_BoolModel):
    enable: bool = False
    sources: TestitemMoleculeEnrichmentSourcesCfg = Field(
        default_factory=lambda: TestitemMoleculeEnrichmentSourcesCfg()
    )
    output: TestitemMoleculeEnrichmentOutputCfg = Field(
        default_factory=lambda: TestitemMoleculeEnrichmentOutputCfg()
    )
    flags: TestitemMoleculeEnrichmentFlagsCfg = Field(
        default_factory=lambda: TestitemMoleculeEnrichmentFlagsCfg()
    )
    logging: TestitemMoleculeEnrichmentLoggingCfg = Field(
        default_factory=lambda: TestitemMoleculeEnrichmentLoggingCfg()
    )

    @field_validator("enable", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class DocumentPubmedCfg(_BaseModel):
    column: str = "PMID"
    sleep: float = Field(5.0, ge=0)
    workers: int = Field(1, ge=1)
    batch_size: int = Field(100, ge=1)
    limit: int | None = Field(default=None, ge=0)


class DocumentChemblCfg(_BaseModel):
    column: str = "document_chembl_id"
    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)


class DocumentAllCfg(_BaseModel):
    column: str = "document_chembl_id"
    chunk_size: int = Field(5, ge=1)
    sleep: float = Field(5.0, ge=0)
    workers: int = Field(1, ge=1)
    batch_size: int = Field(50, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)


class DocumentCfg(_BaseModel):
    pubmed: DocumentPubmedCfg = Field(default_factory=lambda: DocumentPubmedCfg())
    chembl: DocumentChemblCfg = Field(default_factory=lambda: DocumentChemblCfg())
    all: DocumentAllCfg = Field(default_factory=lambda: DocumentAllCfg())


class TargetUniprotCfg(_BaseModel):
    column: str = "uniprot_id"
    data_dir: Path = Path("dictionary/_target/_uniprot")
    limit: int | None = Field(default=None, ge=0)


class TargetChemblCfg(_BaseModel):
    """Defaults for fetching ChEMBL targets.

    The ``column`` aligns with the exported ``target_chembl_id`` identifier
    used throughout the pipelines and CLI defaults.
    """

    column: str = "target_chembl_id"

    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)
    limit: int | None = Field(default=None, ge=0)


class TargetIupharCfg(_BaseModel):
    target_csv: Path = Path("dictionary/_target/_IUPHAR/_IUPHAR_target.csv")
    family_csv: Path = Path("dictionary/_target/_IUPHAR/_IUPHAR_family.csv")
    limit: int | None = Field(default=None, ge=0)


class TargetAllCfg(_BaseModel):
    data_dir: Path = Path("dictionary/_target/_uniprot")
    target_csv: Path = Path("dictionary/_target/_IUPHAR/_IUPHAR_target.csv")
    family_csv: Path = Path("dictionary/_target/_IUPHAR/_IUPHAR_family.csv")
    chunk_size: int = Field(5, ge=1)
    timeout: float = Field(30.0, ge=0)
    uniprot_column: str = "uniprot_id"
    chembl_out: Path | None = None
    uniprot_out: Path | None = None
    iuphar_out: Path | None = None
    limit: int | None = Field(default=None, ge=0)


class TargetCfg(_BaseModel):
    uniprot: TargetUniprotCfg = Field(default_factory=lambda: TargetUniprotCfg())
    chembl: TargetChemblCfg = Field(default_factory=lambda: TargetChemblCfg())
    iuphar: TargetIupharCfg = Field(default_factory=lambda: TargetIupharCfg())
    all: TargetAllCfg = Field(default_factory=lambda: TargetAllCfg())


class ChemblPipelinesCfg(_BaseModel):
    activity: ActivityCfg = Field(default_factory=lambda: ActivityCfg())
    assay: AssayCfg = Field(default_factory=lambda: AssayCfg())
    testitem: TestitemCfg = Field(default_factory=lambda: TestitemCfg())
    document: DocumentCfg = Field(default_factory=lambda: DocumentCfg())
    target: TargetCfg = Field(default_factory=lambda: TargetCfg())


class ChemblSourceCfg(_BaseModel):
    api: ApiCfg = Field(default_factory=lambda: ApiCfg())
    cache: ChemblCacheCfg = Field(default_factory=lambda: ChemblCacheCfg())
    molecule_catalog: MoleculeCatalogCfg = Field(
        default_factory=lambda: MoleculeCatalogCfg()
    )
    pipelines: ChemblPipelinesCfg = Field(default_factory=lambda: ChemblPipelinesCfg())


class UniprotSourceCfg(_BaseModel):
    api: UniprotCfg = Field(default_factory=lambda: UniprotCfg())
    mapping: UniprotMappingCfg = Field(default_factory=lambda: UniprotMappingCfg())


class SourcesCfg(_BaseModel):
    chembl: ChemblSourceCfg = Field(default_factory=lambda: ChemblSourceCfg())
    openalex: OpenAlexCfg = Field(default_factory=lambda: OpenAlexCfg())
    crossref: CrossRefCfg = Field(default_factory=lambda: CrossRefCfg())
    uniprot: UniprotSourceCfg = Field(default_factory=lambda: UniprotSourceCfg())
    iuphar: IupharCfg = Field(default_factory=lambda: IupharCfg())
    pubchem: PubChemCfg = Field(default_factory=lambda: PubChemCfg())
    pubmed: PubMedCfg = Field(default_factory=lambda: PubMedCfg())
    semantic_scholar: SemanticScholarCfg = Field(
        default_factory=lambda: SemanticScholarCfg()
    )


class LocalCfg(_BaseModel):
    resources: ResourcesCfg = Field(default_factory=lambda: ResourcesCfg())
    io: IoCfg = Field(default_factory=lambda: IoCfg())
    init: InitCfg = Field(default_factory=lambda: InitCfg())


class SystemCfg(_BaseModel):
    log: LogCfg = Field(default_factory=lambda: LogCfg())
    rate: RateCfg = Field(default_factory=lambda: RateCfg())
    retry: RetryCfg = Field(default_factory=lambda: RetryCfg())
    doc_type: DocTypeCfg = Field(default_factory=lambda: DocTypeCfg())
    doc_quality: DocQualityCfg = Field(default_factory=lambda: DocQualityCfg())


class Config(_BaseModel):
    sources: SourcesCfg = Field(default_factory=lambda: SourcesCfg())
    local: LocalCfg = Field(default_factory=lambda: LocalCfg())
    system: SystemCfg = Field(default_factory=lambda: SystemCfg())
    activity_bounds: ActivityBoundsCfg = Field(
        default_factory=lambda: ActivityBoundsCfg()
    )
    activity_enrichment: ActivityEnrichmentCfg = Field(
        default_factory=lambda: ActivityEnrichmentCfg()
    )
    testitem_molecule_enrichment: TestitemMoleculeEnrichmentCfg = Field(
        default_factory=lambda: TestitemMoleculeEnrichmentCfg()
    )

    # -- Compatibility accessors -------------------------------------------------
    #
    # The project historically exposed flat attributes such as ``cfg.api`` and
    # ``cfg.io``.  The configuration structure is now grouped by responsibility,
    # but the properties below keep the public interface stable while we update
    # callers incrementally.

    @property
    def api(self) -> ApiCfg:
        return self.sources.chembl.api

    @property
    def chembl(self) -> ChemblCacheCfg:
        return self.sources.chembl.cache

    @property
    def molecule_catalog(self) -> MoleculeCatalogCfg:
        return self.sources.chembl.molecule_catalog

    @property
    def openalex(self) -> OpenAlexCfg:
        return self.sources.openalex

    @property
    def crossref(self) -> CrossRefCfg:
        return self.sources.crossref

    @property
    def uniprot(self) -> UniprotCfg:
        return self.sources.uniprot.api

    @property
    def uniprot_mapping(self) -> UniprotMappingCfg:
        return self.sources.uniprot.mapping

    @property
    def iuphar(self) -> IupharCfg:
        return self.sources.iuphar

    @property
    def pubchem(self) -> PubChemCfg:
        return self.sources.pubchem

    @property
    def pubmed(self) -> PubMedCfg:
        return self.sources.pubmed

    @property
    def semantic_scholar(self) -> SemanticScholarCfg:
        return self.sources.semantic_scholar

    @property
    def doc_type(self) -> DocTypeCfg:
        return self.system.doc_type

    @property
    def doc_quality(self) -> DocQualityCfg:
        return self.system.doc_quality

    @property
    def resources(self) -> ResourcesCfg:
        return self.local.resources

    @property
    def io(self) -> IoCfg:
        return self.local.io

    @property
    def log(self) -> LogCfg:
        return self.system.log

    @property
    def init(self) -> InitCfg:
        return self.local.init

    @property
    def rate(self) -> RateCfg:
        return self.system.rate

    @property
    def retry(self) -> RetryCfg:
        return self.system.retry

    @property
    def activity(self) -> ActivityCfg:
        return self.sources.chembl.pipelines.activity

    @property
    def assay(self) -> AssayCfg:
        return self.sources.chembl.pipelines.assay

    @property
    def testitem(self) -> TestitemCfg:
        return self.sources.chembl.pipelines.testitem

    @property
    def document(self) -> DocumentCfg:
        return self.sources.chembl.pipelines.document

    @property
    def target(self) -> TargetCfg:
        return self.sources.chembl.pipelines.target


# Dotted paths for configuration fields that are defined as ``Path`` types.
_CONFIG_PATH_FIELDS: set[tuple[str, ...]] = _collect_path_field_paths(Config)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------


def session_with_retry(api: ApiCfg, retry: RetryCfg) -> Session:
    """Return an HTTP session configured for retries and user agent.

    The returned session retries failed requests for *all* HTTP methods,
    including ``POST``. Requests are attempted once plus
    ``retry.max_attempts - 1`` automatic retries for retryable statuses. The
    session avoids raising exceptions on HTTP error status codes, allowing
    callers to handle responses manually.

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
        total=max(0, retry.max_attempts - 1),
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


def _session_with_mailto_header(api: ApiCfg, retry: RetryCfg, mailto: str) -> Session:
    """Return a session configured with the ``mailto`` contact header."""

    session = session_with_retry(api, retry)
    session.headers["mailto"] = mailto
    return session


def openalex_session(api: ApiCfg, retry: RetryCfg, cfg: OpenAlexCfg) -> Session:
    """Return a session configured for OpenAlex requests."""

    return _session_with_mailto_header(api, retry, cfg.mailto)


def crossref_session(api: ApiCfg, retry: RetryCfg, cfg: CrossRefCfg) -> Session:
    """Return a session configured for CrossRef requests."""

    return _session_with_mailto_header(api, retry, cfg.mailto)


def _set_by_path(data: dict[str, Any], path: list[str], value: Any) -> None:
    cur: dict[str, Any] = data
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[path[-1]] = value


def _apply_env_overrides(data: dict[str, Any]) -> dict[tuple[str, ...], str]:
    prefix = "CHEMBL_DA"
    overrides: dict[tuple[str, ...], str] = {}
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
        value = _parse_env_value(key, env_val)
        _set_by_path(data, parts, value)
        overrides[tuple(parts)] = key
    return overrides


def _parse_env_value(env_key: str, raw_value: str) -> Any:
    """Normalize *raw_value* from environment variable *env_key*.

    Empty strings are returned unchanged to avoid coercing intentional blanks to
    ``None``. Other values are parsed using :func:`yaml.safe_load` so that
    numbers, booleans, lists, and mappings are deserialized before validation.
    """

    if raw_value == "":
        return ""
    if raw_value and raw_value.strip() == "":
        return raw_value
    try:
        value = yaml.safe_load(raw_value)
    except yaml.YAMLError as exc:
        logger.debug(
            "treating %s as plain string due to YAML parse error: %s",
            env_key,
            exc,
        )
        return raw_value
    return _coerce_integral_numbers(value)


def _normalize_env_errors(
    errors: Sequence[ErrorDetails], overrides: Mapping[tuple[str, ...], str]
) -> tuple[list[str], int]:
    """Return formatted messages for validation *errors* caused by overrides."""

    messages: list[str] = []
    handled = 0
    for error in errors:
        loc = error.get("loc", ())
        if not isinstance(loc, Sequence):
            continue
        str_path = [str(part).lower() for part in loc if isinstance(part, str)]
        env_key: str | None = None
        for index in range(len(str_path), 0, -1):
            candidate = tuple(str_path[:index])
            match = overrides.get(candidate)
            if match:
                env_key = match
                break
        if not env_key:
            continue
        messages.append(_format_env_error(env_key, error))
        handled += 1
    return messages, handled


def _format_env_error(env_key: str, error: Mapping[str, Any]) -> str:
    """Return a human readable error string for *env_key* based on *error*."""

    message = _format_env_error_message(error)
    loc = error.get("loc")
    location = _format_error_location(loc) if isinstance(loc, Sequence) else ""
    if location:
        return f"{env_key} ({location}) {message}".strip()
    return f"{env_key} {message}".strip()


def _format_env_error_message(error: Mapping[str, Any]) -> str:
    """Convert a Pydantic error dictionary into a concise human message."""

    error_type = error.get("type")
    ctx: dict[str, Any] = error.get("ctx") or {}
    if error_type == "greater_than_equal" and "ge" in ctx:
        return f"must be ≥{ctx['ge']}"
    if error_type == "greater_than" and "gt" in ctx:
        return f"must be >{ctx['gt']}"
    if error_type == "less_than_equal" and "le" in ctx:
        return f"must be ≤{ctx['le']}"
    if error_type == "less_than" and "lt" in ctx:
        return f"must be <{ctx['lt']}"
    if error_type == "int_parsing":
        return "must be an integer"
    if error_type == "float_parsing":
        return "must be a number"
    if error_type == "bool_parsing":
        return "must be a boolean"
    if error_type in {"string_type", "string_parsing"}:
        return "must be a string"
    if error_type == "list_type":
        return "must be a list"
    if error_type == "dict_type":
        return "must be a mapping"
    if error_type in {"enum", "literal_error"} and "expected" in ctx:
        expected = ", ".join(map(str, ctx["expected"]))
        return f"must be one of {expected}"
    raw_message = error.get("msg")
    if isinstance(raw_message, str):
        message = raw_message
    elif raw_message is None:
        message = "is invalid"
    else:
        message = str(raw_message)
    if message.startswith("Input should be "):
        return "must be " + message[len("Input should be ") :]
    return message


def _format_error_location(loc: Sequence[Any] | None) -> str:
    if not loc:
        return ""
    parts: list[str] = []
    for part in loc:
        if isinstance(part, int):
            if parts:
                parts[-1] = f"{parts[-1]}[{part}]"
            else:
                parts.append(f"[{part}]")
        else:
            parts.append(str(part))
    return ".".join(parts)


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


_OPTIONAL_UNKNOWN_KEYS: frozenset[str] = frozenset(
    {
        "local.io.csv_fallback_separators",
    }
)


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
    path: str | Path | None = None,
    cli_overrides: dict[str, Any] | None = None,
    *,
    base_path: Path | str | None = None,
    strict: bool = True,
) -> Config:
    """Load configuration from *path* applying overrides.

    Parameters
    ----------
    base_path:
        Optional root directory used to resolve ``$CHEMBL_DA_BASE_PATH``
        placeholders.
    """

    try:
        data, resolved_path = load_yaml_config(path)
    except ConfigLoaderError as exc:
        raise ConfigError(str(exc)) from exc

    local_path = resolved_path.with_name(
        f"{resolved_path.stem}.local{resolved_path.suffix}"
    )
    if local_path != resolved_path and local_path.exists():
        try:
            local_data, _ = load_yaml_config(local_path)
        except ConfigLoaderError as exc:
            raise ConfigError(str(exc)) from exc
        _merge_mapping(data, local_data)

    data = _coerce_integral_numbers(data)

    # Guard against accidentally passing the JSON schema instead of a runtime
    # configuration file. The schema contains the ``$defs`` key at the top
    # level which is not present in actual application settings.
    if "$defs" in data:
        raise ConfigError(
            f"{resolved_path} appears to be a configuration schema; "
            "provide an application config file such as config/config.yaml."
        )

    env_overrides = _apply_env_overrides(data)
    _upgrade_legacy_config(data)

    if cli_overrides:
        for key, val in cli_overrides.items():
            _set_by_path(data, key.split("."), val)

    placeholder_base = _resolve_placeholder_base_path(base_path)
    data = _expand_config_placeholders(data, base_path=placeholder_base)

    base_dir = resolved_path.parent.resolve()
    _absolutise_config_paths(data, base_dir)

    unknown = _collect_unknown_keys(data, Config)
    ignored_unknown = [
        key for key in unknown if key in _OPTIONAL_UNKNOWN_KEYS
    ]
    unknown = [key for key in unknown if key not in _OPTIONAL_UNKNOWN_KEYS]
    if ignored_unknown:
        logger.warning(
            "config_unknown_ignored",
            keys=", ".join(sorted(ignored_unknown)),
            hint="Upgrade the application to use these options.",
        )
    if unknown:
        msg = (
            "Unknown configuration key(s) in "
            f"{resolved_path}: {', '.join(sorted(unknown))}"
        )
        if strict:
            raise ValueError(msg)
        logger.warning(msg)

    try:
        cfg = Config.model_validate(data)
    except ValidationError as exc:
        env_messages, handled = _normalize_env_errors(exc.errors(), env_overrides)
        if env_messages:
            message = "; ".join(env_messages)
            if handled < len(exc.errors()):
                message = f"{message}; additional validation errors: {exc}"
            raise ConfigError(message) from exc
        raise

    if not cfg.io.exist_ok:
        for p in (cfg.io.output_dir, cfg.io.cache_dir):
            if not p.exists():
                raise FileNotFoundError(p)

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


def _merge_mapping(dest: dict[str, Any], src: dict[str, Any]) -> None:
    """Recursively merge mapping *src* into *dest*."""

    for key, value in src.items():
        if key in dest and isinstance(dest[key], dict) and isinstance(value, dict):
            _merge_mapping(dest[key], value)
        else:
            dest[key] = value


def _upgrade_legacy_config(data: dict[str, Any]) -> None:
    """Translate legacy flat config keys into the structured layout."""

    sources = data.setdefault("sources", {})
    chembl = sources.setdefault("chembl", {})
    pipelines = chembl.setdefault("pipelines", {})
    local = data.setdefault("local", {})
    system_cfg = data.setdefault("system", {})

    if "api" in data:
        chembl.setdefault("api", {})
        _merge_mapping(chembl["api"], data.pop("api"))
    if "chembl" in data:
        chembl.setdefault("cache", {})
        _merge_mapping(chembl["cache"], data.pop("chembl"))

    for section in (
        "openalex",
        "crossref",
        "iuphar",
        "pubchem",
        "pubmed",
        "semantic_scholar",
    ):
        if section in data:
            sources.setdefault(section, {})
            _merge_mapping(sources[section], data.pop(section))

    if "uniprot" in data:
        sources.setdefault("uniprot", {})
        api_cfg = sources["uniprot"].setdefault("api", {})
        _merge_mapping(api_cfg, data.pop("uniprot"))
    if "uniprot_mapping" in data:
        sources.setdefault("uniprot", {})
        mapping_cfg = sources["uniprot"].setdefault("mapping", {})
        _merge_mapping(mapping_cfg, data.pop("uniprot_mapping"))

    for section in ("activity", "assay", "testitem", "document", "target"):
        if section in data:
            pipelines.setdefault(section, {})
            _merge_mapping(pipelines[section], data.pop(section))

    if "resources" in data:
        local.setdefault("resources", {})
        _merge_mapping(local["resources"], data.pop("resources"))
    if "io" in data:
        local.setdefault("io", {})
        _merge_mapping(local["io"], data.pop("io"))
    if "init" in data:
        local.setdefault("init", {})
        _merge_mapping(local["init"], data.pop("init"))

    if "log" in data:
        system_cfg.setdefault("log", {})
        _merge_mapping(system_cfg["log"], data.pop("log"))
    if "rate" in data:
        system_cfg.setdefault("rate", {})
        _merge_mapping(system_cfg["rate"], data.pop("rate"))
    if "retry" in data:
        system_cfg.setdefault("retry", {})
        _merge_mapping(system_cfg["retry"], data.pop("retry"))
    if "doc_type" in data:
        system_cfg.setdefault("doc_type", {})
        _merge_mapping(system_cfg["doc_type"], data.pop("doc_type"))


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
    "CHEMBL_DA_BASE": ["sources", "chembl", "api", "chembl_base"],
    "CHEMBL_DA_BURST": ["sources", "chembl", "api", "burst"],
    "CHEMBL_DA_CACHE_DIR": ["local", "io", "cache_dir"],
    "CHEMBL_DA__IO__CACHE_DIR": ["local", "io", "cache_dir"],
    "CHEMBL_DA__IO__EXIST_OK": ["local", "io", "exist_ok"],
    "CHEMBL_DA_CACHE_MAXSIZE": ["sources", "chembl", "cache", "cache_maxsize"],
    "CHEMBL_DA_CACHE_TTL": ["sources", "chembl", "cache", "cache_ttl"],
    "CHEMBL_DA_MOLECULE_CATALOG_CACHE": [
        "sources",
        "chembl",
        "molecule_catalog",
        "cache_path",
    ],
    "CHEMBL_DA_DICT_DIR": ["local", "resources", "dictionary_dir"],
    "CHEMBL_DA_GLOBAL_BURST": ["system", "rate", "global_burst"],
    "CHEMBL_DA_GLOBAL_RPS": ["system", "rate", "global_rps"],
    "CHEMBL_DA_IUPHAR_FAMILY_CSV": ["local", "resources", "iuphar_family_csv"],
    "CHEMBL_DA_IUPHAR_TARGET_CSV": ["local", "resources", "iuphar_target_csv"],
    "CHEMBL_DA_LIMITER_CACHE_MAXSIZE": ["system", "rate", "limiter_cache_maxsize"],
    "CHEMBL_DA_LIMITER_CACHE_TTL": ["system", "rate", "limiter_cache_ttl"],
    "CHEMBL_DA_OUTDIR": ["local", "io", "output_dir"],
    "CHEMBL_DA_RPS": ["sources", "chembl", "api", "rps"],
    "CHEMBL_DA_PUBMED_RPS": ["sources", "pubmed", "rps"],
    "CHEMBL_DA_PUBMED_BURST": ["sources", "pubmed", "burst"],
    "CHEMBL_DA_SEMANTIC_SCHOLAR_RPS": ["sources", "semantic_scholar", "rps"],
    "CHEMBL_DA_SEMANTIC_SCHOLAR_BURST": [
        "sources",
        "semantic_scholar",
        "burst",
    ],
    "CHEMBL_DA_TARGETS_TYPE_CSV": ["local", "resources", "targets_type_csv"],
    "CHEMBL_DA_TIMEOUT_CONNECT": ["sources", "chembl", "api", "timeout_connect"],
    "CHEMBL_DA_TIMEOUT_READ": ["sources", "chembl", "api", "timeout_read"],
    "CHEMBL_DA_UNIPROT_DATA_DIR": ["local", "resources", "uniprot_data_dir"],
    "CHEMBL_DA_OPENALEX_MAILTO": ["sources", "openalex", "mailto"],
    "CHEMBL_DA_CROSSREF_MAILTO": ["sources", "crossref", "mailto"],
    "CHEMBL_DA_OPENALEX_TIMEOUT_CONNECT": [
        "sources",
        "openalex",
        "timeout_connect",
    ],
    "CHEMBL_DA_OPENALEX_TIMEOUT_READ": [
        "sources",
        "openalex",
        "timeout_read",
    ],
    "CHEMBL_DA_CROSSREF_TIMEOUT_CONNECT": [
        "sources",
        "crossref",
        "timeout_connect",
    ],
    "CHEMBL_DA_CROSSREF_TIMEOUT_READ": [
        "sources",
        "crossref",
        "timeout_read",
    ],
    "CHEMBL_DA_OPENALEX_BASE": ["sources", "openalex", "base"],
    "CHEMBL_DA_CROSSREF_BASE": ["sources", "crossref", "base"],
    "CHEMBL_DA_UNIPROT_BASE": ["sources", "uniprot", "api", "base"],
    "CHEMBL_DA_IUPHAR_BASE": ["sources", "iuphar", "base"],
    "CHEMBL_DA_PUBCHEM_BASE": ["sources", "pubchem", "base"],
    "CHEMBL_DA_PUBCHEM_USER_AGENT": ["sources", "pubchem", "user_agent"],
    "CHEMBL_DA_LOG_LEVEL": ["system", "log", "level"],
    "CHEMBL_DA_RETRY_MAX_ATTEMPTS": ["system", "retry", "max_attempts"],
    "CHEMBL_DA_RETRY_BACKOFF_FACTOR": ["system", "retry", "backoff_factor"],
    "CHEMBL_DA_RETRY_BACKOFF_CAP": ["system", "retry", "backoff_cap"],
}

_ALIAS_MAP: dict[str, list[str]] = {
    **build_alias_map(Config),
    **_ALIAS_OVERRIDES,
}


__all__ = [
    "ApiCfg",
    "ChemblCacheCfg",
    "MoleculeCatalogCfg",
    "OpenAlexCfg",
    "CrossRefCfg",
    "UniprotCfg",
    "UniprotMappingCfg",
    "IupharCfg",
    "PubChemCfg",
    "PubMedCfg",
    "SemanticScholarCfg",
    "DocTypeCfg",
    "DocQualityCfg",
    "ActivityCfg",
    "ActivityActionTypeCfg",
    "ActivityEnrichmentCfg",
    "ActivityPropertiesCfg",
    "AssayCfg",
    "TestitemCfg",
    "TestitemBatchRetryCfg",
    "TestitemMoleculeEnrichmentCfg",
    "TestitemMoleculeEnrichmentFlagsCfg",
    "TestitemMoleculeEnrichmentLoggingCfg",
    "TestitemMoleculeEnrichmentOutputCfg",
    "TestitemMoleculeEnrichmentSourcesCfg",
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
    "InitCfg",
    "RateCfg",
    "RetryCfg",
    "session_with_retry",
    "openalex_session",
    "crossref_session",
    "LogCfg",
    "ChemblSourceCfg",
    "SourcesCfg",
    "LocalCfg",
    "SystemCfg",
    "Config",
    "ConfigError",
    "ensure_dirs",
    "load_config",
    "print_config",
    "build_alias_map",
]
