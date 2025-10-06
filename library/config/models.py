"""Pydantic configuration schemas for the ChEMBL data acquisition tool."""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence
from urllib.parse import urlparse

from pydantic import BaseModel, ConfigDict, Field, StrictInt, field_validator

from config import DICTIONARY_DIR

from ..document_defaults import ALL_DEFAULTS, CHEMBL_DEFAULTS, PUBMED_DEFAULTS














def _resolve_default_path(value: Path | str) -> Path:
    candidate = Path(value).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    return (Path.cwd() / candidate).resolve()


def _default_base_path() -> Path:
    """Return the default data directory used for placeholder expansion."""

    env_override = os.environ.get("CHEMBL_DA_BASE_PATH")
    if env_override:
        return _resolve_default_path(env_override)
    return (Path.home() / ".local" / "share" / "chembl-da").resolve()


def _default_cache_home() -> Path:
    """Return the default cache directory for local artefacts."""

    return (Path.home() / ".cache" / "chembl-da").resolve()






_UNSET = object()


@dataclass(slots=True)
class ConfigSource:
    """Origin of a configuration value used in :class:`ConfigMetadata`."""

    source: str
    detail: str | None = None


@dataclass(slots=True)
class ConfigMetadata:
    """Metadata describing configuration values and their provenance."""

    snapshot: dict[str, Any]
    sources: dict[tuple[str, ...], ConfigSource]
    cli_paths: dict[str, tuple[str, ...]] = field(default_factory=dict)

    def get(self, path: Sequence[str] | str) -> dict[str, Any]:
        """Return a copy of the metadata entry located at ``path``."""

        tuple_path = (
            tuple(str(part) for part in path)
            if isinstance(path, Sequence) and not isinstance(path, str)
            else tuple(str(part) for part in str(path).split(".") if str(part))
        )
        current: Any = self.snapshot
        for key in tuple_path:
            if not isinstance(current, dict) or key not in current:
                return {"value": None, "source": "unknown"}
            current = current[key]
        if isinstance(current, dict) and "value" in current and "source" in current:
            result = {"value": current["value"], "source": current["source"]}
            if "detail" in current and current["detail"] is not None:
                result["detail"] = current["detail"]
            return result
        return {"value": current, "source": "unknown"}

    def cli_entry(self, argument: str) -> dict[str, Any] | None:
        """Return metadata for the CLI destination ``argument`` if available."""

        path = self.cli_paths.get(argument)
        if path is None:
            return None
        return self.get(path)

    def option(
        self,
        *,
        argument: str | None = None,
        path: Sequence[str] | str | None = None,
        value: Any = _UNSET,
        default_source: str = "unknown",
        default_detail: str | None = None,
    ) -> dict[str, Any]:
        """Return metadata for ``argument`` or ``path`` overriding ``value``."""

        entry: dict[str, Any] | None = None
        if argument is not None:
            entry = self.cli_entry(argument)
        if entry is None and path is not None:
            entry = self.get(path)
        if entry is None:
            entry = {"value": None, "source": default_source}
            if default_detail is not None:
                entry["detail"] = default_detail
        else:
            entry = dict(entry)
        if value is not _UNSET:
            entry["value"] = value
        return entry






def _dictionary_resource(*parts: str) -> Path:
    """Return a filesystem path for a bundled dictionary resource."""

    return DICTIONARY_DIR.joinpath(*parts)


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
    timeout_read: float = Field(60.0, ge=1)
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
    hierarchy_lookup_path: Path | None = (
        DICTIONARY_DIR / "_testitem" / "molecule_hierarchy.csv"
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
    max_attempts: int = Field(4, ge=1)
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
    timeout: float = Field(90.0, gt=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    dry_run: bool = False

    @field_validator("dry_run", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class AssayCfg(_BaseModel):
    column: str = "assay_chembl_id"
    batch_size: int = Field(25, ge=1)
    timeout: float = Field(60.0, gt=0)
    limit: int | None = Field(default=None, ge=0)


class CelllineCfg(_BaseModel):
    column: str = "cell_chembl_id"
    batch_size: int = Field(20, ge=1)
    timeout: float = Field(30.0, gt=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)


class TissueCfg(_BaseModel):
    column: str = "tissue_chembl_id"
    batch_size: int = Field(20, ge=1)
    timeout: float = Field(30.0, gt=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)


class TestitemBatchRetryCfg(_BoolModel):
    enable: bool = True
    shrink_factor: float = Field(0.5, gt=0, lt=1)
    min_size: int = Field(1, ge=1)

    @field_validator("enable", mode="before")
    @classmethod
    def _bools(cls, v: Any) -> bool:
        return cls._parse_bool(v)


class TestitemCfg(_BaseModel):
    column: str = "molecule_chembl_id"
    batch_size: int = Field(250, ge=1, le=1000)
    timeout: float = Field(90.0, gt=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    fields: tuple[str, ...] = Field(default=TESTITEM_FIELD_DEFAULTS)
    request_limit: int = Field(1000, ge=1, le=1000)
    retries: int | None = Field(default=None, ge=0)
    backoff_factor: float | None = Field(default=None, ge=0)
    batch_retry: TestitemBatchRetryCfg = Field(
        default_factory=lambda: TestitemBatchRetryCfg()
    )
    parent_watchdog_idle_minutes: float = Field(0.0, ge=0)
    execution_budget_minutes: float | None = Field(default=None, ge=0)

    @field_validator("fields", mode="before")
    @classmethod
    def _coerce_fields(cls, value: Any) -> tuple[str, ...]:
        if value is None:
            return TESTITEM_FIELD_DEFAULTS
        if isinstance(value, (list, tuple)):
            return tuple(str(item) for item in value)
        raise TypeError("testitem.fields must be a list or tuple of strings")


class TestitemMoleculeEnrichmentSourcesCfg(_BaseModel):
    molecule_catalog_path: Path = (
        DICTIONARY_DIR / "_testitem" / "molecule_catalog.csv"
    )
    molecule_hierarchy_path: Path = (
        DICTIONARY_DIR / "_testitem" / "molecule_hierarchy.csv"
    )


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
    column: str = PUBMED_DEFAULTS.column
    sleep: float = Field(PUBMED_DEFAULTS.sleep, ge=0)
    workers: int = Field(PUBMED_DEFAULTS.workers, ge=1)
    batch_size: int = Field(PUBMED_DEFAULTS.batch_size, ge=1)
    limit: int | None = Field(default=None, ge=0)


class DocumentChemblCfg(_BaseModel):
    column: str = CHEMBL_DEFAULTS.column
    chunk_size: int = Field(CHEMBL_DEFAULTS.chunk_size, ge=1)
    timeout: float = Field(CHEMBL_DEFAULTS.timeout, gt=0)
    limit: int | None = Field(default=None, ge=0)


class DocumentAllCfg(_BaseModel):
    column: str = ALL_DEFAULTS.column
    chunk_size: int = Field(ALL_DEFAULTS.chunk_size, ge=1)
    sleep: float = Field(ALL_DEFAULTS.sleep, ge=0)
    workers: int = Field(ALL_DEFAULTS.workers, ge=1)
    batch_size: int = Field(ALL_DEFAULTS.batch_size, ge=1)
    timeout: float = Field(ALL_DEFAULTS.timeout, gt=0)
    limit: int | None = Field(default=None, ge=0)


class DocumentCfg(_BaseModel):
    pubmed: DocumentPubmedCfg = Field(default_factory=lambda: DocumentPubmedCfg())
    chembl: DocumentChemblCfg = Field(default_factory=lambda: DocumentChemblCfg())
    all: DocumentAllCfg = Field(default_factory=lambda: DocumentAllCfg())


class TargetUniprotCfg(_BaseModel):
    column: str = "uniprot_id"
    data_dir: Path = DICTIONARY_DIR / "_target" / "_uniprot"
    limit: int | None = Field(default=None, ge=0)
    chunk_size: int = Field(100, ge=1)
    timeout: float = Field(30.0, gt=0)
    offset: int = Field(0, ge=0)


class TargetChemblBatchRetryCfg(_BoolModel):
    enable: bool = True
    shrink_factor: float = Field(0.5, gt=0, lt=1)
    min_size: int = Field(1, ge=1)
    single_timeout_retries: int = Field(2, ge=0)
    single_timeout_delay: float = Field(0.0, ge=0.0)

    @field_validator("enable", mode="before")
    @classmethod
    def _bools(cls, value: Any) -> bool:
        return cls._parse_bool(value)


class TargetChemblCfg(_BaseModel):
    """Defaults for fetching ChEMBL targets.

    The ``column`` aligns with the exported ``target_chembl_id`` identifier
    used throughout the pipelines and CLI defaults.
    """

    column: str = "target_chembl_id"

    chunk_size: int = Field(3, ge=1)
    timeout: float = Field(90.0, gt=0)

    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    batch_retry: TargetChemblBatchRetryCfg = Field(
        default_factory=lambda: TargetChemblBatchRetryCfg()
    )



class TargetIupharCfg(_BaseModel):
    column: str = "uniprot_id"
    chunk_size: int = Field(100, ge=1)
    timeout: float = Field(30.0, gt=0)
    limit: int | None = Field(default=None, ge=0)
    offset: int = Field(0, ge=0)
    target_csv: Path = (
        DICTIONARY_DIR / "_target" / "_IUPHAR" / "_IUPHAR_target.csv"
    )
    family_csv: Path = (
        DICTIONARY_DIR / "_target" / "_IUPHAR" / "_IUPHAR_family.csv"
    )


class TargetAllCfg(_BaseModel):
    """Defaults for the combined target pipeline."""

    column: str = "target_chembl_id"
    data_dir: Path = DICTIONARY_DIR / "_target" / "_uniprot"
    target_csv: Path = (
        DICTIONARY_DIR / "_target" / "_IUPHAR" / "_IUPHAR_target.csv"
    )
    family_csv: Path = (
        DICTIONARY_DIR / "_target" / "_IUPHAR" / "_IUPHAR_family.csv"
    )
    chunk_size: int = Field(3, ge=1)
    timeout: float = Field(90.0, gt=0)
    uniprot_column: str = "uniprot_id"
    chembl_out: Path | None = None
    uniprot_out: Path | None = None
    iuphar_out: Path | None = None
    limit: int | None = Field(default=None, ge=1)
    offset: int = Field(0, ge=0)


class TargetCfg(_BaseModel):
    uniprot: TargetUniprotCfg = Field(default_factory=lambda: TargetUniprotCfg())
    chembl: TargetChemblCfg = Field(default_factory=lambda: TargetChemblCfg())
    iuphar: TargetIupharCfg = Field(default_factory=lambda: TargetIupharCfg())
    all: TargetAllCfg = Field(default_factory=lambda: TargetAllCfg())


class ChemblPipelinesCfg(_BaseModel):
    activity: ActivityCfg = Field(default_factory=lambda: ActivityCfg())
    assay: AssayCfg = Field(default_factory=lambda: AssayCfg())
    cellline: CelllineCfg = Field(default_factory=lambda: CelllineCfg())
    tissue: TissueCfg = Field(default_factory=lambda: TissueCfg())
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
    def cellline(self) -> CelllineCfg:
        return self.sources.chembl.pipelines.cellline

    @property
    def tissue(self) -> TissueCfg:
        return self.sources.chembl.pipelines.tissue

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


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

















































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
    "ActivityBoundsCfg",
    "ActivityCfg",
    "ActivityActionTypeCfg",
    "ActivityEnrichmentCfg",
    "ActivityPropertiesCfg",
    "AssayCfg",
    "TestitemCfg",
    "TestitemBatchRetryCfg",
    "TESTITEM_FIELD_DEFAULTS",
    "TestitemMoleculeEnrichmentCfg",
    "TestitemMoleculeEnrichmentFlagsCfg",
    "TestitemMoleculeEnrichmentLoggingCfg",
    "TestitemMoleculeEnrichmentOutputCfg",
    "TestitemMoleculeEnrichmentSourcesCfg",
    "DocumentPubmedCfg",
    "DocumentChemblCfg",
    "DocumentAllCfg",
    "ConfigMetadata",
    "ConfigSource",
    "DocumentCfg",
    "TargetUniprotCfg",
    "TargetChemblBatchRetryCfg",
    "TargetChemblCfg",
    "TargetIupharCfg",
    "TargetAllCfg",
    "TargetCfg",
    "ResourcesCfg",
    "IoCfg",
    "InitCfg",
    "RateCfg",
    "RetryCfg",
    "LogCfg",
    "ChemblSourceCfg",
    "SourcesCfg",
    "LocalCfg",
    "SystemCfg",
    "Config",
    "ConfigError",
]
