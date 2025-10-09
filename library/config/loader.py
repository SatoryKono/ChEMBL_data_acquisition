"""Helpers for loading and normalising configuration files."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Mapping, Sequence
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping as TypingMapping

import yaml
from pydantic import BaseModel, ValidationError

from ..common.log import logger
from library.resources.dictionaries import (
    DictionaryManifestError,
    list_resource_names,
    list_resources,
)
from config.paths import CONFIG_DIR as _CONFIG_DIR
from config.paths import DEFAULT_CONFIG_PATH as _DEFAULT_CONFIG_PATH
from .runtime import configure_rate_limiters
from .env import _apply_env_overrides
from .env import _expand_config_placeholders
from .env import _normalize_env_errors
from .env import _resolve_placeholder_base_path
from .env import _set_by_path
from .models import Config, ConfigError, ConfigMetadata, ConfigSource, _CONFIG_PATH_FIELDS

__all__ = [
    "CONFIG_DIR",
    "DEFAULT_CONFIG_PATH",
    "DEFAULT_CONFIG_RELATIVE",
    "ConfigLoaderError",
    "ConfigMetadata",
    "ConfigSource",
    "ConfigError",
    "load_config",
    "load_yaml_config",
    "ensure_dirs",
    "print_config",
    "resolve_config_path",
    "_absolutise_path_value",
    "_serialize_paths",
    "_mask_secrets",
]


CONFIG_DIR = _CONFIG_DIR
DEFAULT_CONFIG_PATH = _DEFAULT_CONFIG_PATH
_DEFAULT_CONFIG_NAME = DEFAULT_CONFIG_PATH.name
DEFAULT_CONFIG_RELATIVE = Path("config") / _DEFAULT_CONFIG_NAME
_PROJECT_ROOT = Path(__file__).resolve().parents[2]


class ConfigLoaderError(RuntimeError):
    """Raised when configuration parsing fails."""


def resolve_config_path(path: str | Path | None = None) -> Path:
    """Return an absolute path to the configuration file."""

    if path is None:
        return DEFAULT_CONFIG_PATH

    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate

    for base in (Path.cwd(), CONFIG_DIR.parent):
        base_candidate = (base / candidate).resolve()
        if base_candidate.exists():
            return base_candidate

    if candidate.parent == Path(".") and candidate.name == _DEFAULT_CONFIG_NAME:
        default_candidate = (CONFIG_DIR / candidate.name).resolve()
        if default_candidate.exists():
            return default_candidate

    return (Path.cwd() / candidate).resolve()


def load_yaml_config(path: str | Path | None = None) -> tuple[dict[str, Any], Path]:
    """Return raw configuration data loaded from YAML."""

    cfg_path = resolve_config_path(path)
    try:
        with cfg_path.open("r", encoding="utf8") as handle:
            data = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:  # pragma: no cover - defensive
        raise ConfigLoaderError(f"configuration file not found: {cfg_path}") from exc
    except yaml.YAMLError as exc:  # pragma: no cover - defensive
        raise ConfigLoaderError(
            f"failed to parse YAML configuration at {cfg_path}: {exc}"
        ) from exc

    if not isinstance(data, dict):
        raise ConfigLoaderError("top-level structure in config file must be a mapping")
    return data, cfg_path


def _format_source_detail(path: Path) -> str:
    """Return a stable, human-readable representation of ``path``."""

    candidate = Path(path)
    try:
        resolved = candidate.resolve()
    except OSError:
        resolved = candidate

    if not resolved.is_absolute():
        return resolved.as_posix()

    bases: list[Path] = []
    try:
        cwd = Path.cwd().resolve()
    except OSError:
        cwd = None
    if cwd is not None:
        bases.append(cwd)

    project_root = _PROJECT_ROOT
    if project_root not in bases:
        bases.append(project_root)

    for base in bases:
        try:
            relative = resolved.relative_to(base)
        except ValueError:
            continue
        if not relative.parts:
            return resolved.name
        return relative.as_posix()

    return resolved.as_posix()


def _absolutise_path_value(value: Any, base_dir: Path) -> Any:
    """Return *value* converted to an absolute path relative to *base_dir*."""

    if value is None:
        return value
    if isinstance(value, (str, os.PathLike)):
        candidate = os.fspath(value)
        if candidate in _dictionary_resource_names():
            return value
    base_dir = base_dir.resolve()

    def _resolve(path: Path) -> Path:
        if path.is_absolute():
            return path
        if base_dir.name and path.parts and path.parts[0] == base_dir.name:
            return (base_dir.parent / path).resolve()
        return (base_dir / path).resolve()

    if isinstance(value, str):
        resolved = _resolve(Path(value))
        return str(resolved)
    if isinstance(value, os.PathLike):
        resolved = _resolve(Path(value))
        return resolved
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


def _iter_leaf_items(
    data: Any, prefix: tuple[str, ...] = ()
) -> Iterator[tuple[tuple[str, ...], Any]]:
    """Yield ``(path, value)`` pairs for leaf nodes in ``data``."""

    if isinstance(data, dict):
        for key, value in data.items():
            key_str = str(key)
            yield from _iter_leaf_items(value, (*prefix, key_str))
        return
    yield prefix, data


def _build_snapshot(
    data: Mapping[str, Any],
    sources: Mapping[tuple[str, ...], ConfigSource],
    prefix: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Return nested mapping annotating ``data`` with source metadata."""

    snapshot: dict[str, Any] = {}
    for key, value in data.items():
        key_str = str(key)
        path = (*prefix, key_str)
        if isinstance(value, Mapping):
            snapshot[key_str] = _build_snapshot(value, sources, path)
            continue
        source = sources.get(path, ConfigSource("default", None))
        entry: dict[str, Any] = {"value": value, "source": source.source}
        if source.detail is not None:
            entry["detail"] = source.detail
        snapshot[key_str] = entry
    return snapshot


_OPTIONAL_UNKNOWN_KEYS: frozenset[str] = frozenset(
    {
        "local.io.csv_fallback_separators",
    }
)


@lru_cache(maxsize=1)
def _dictionary_resource_names() -> frozenset[str]:
    """Return the set of known dictionary resource identifiers."""

    try:
        names = list_resource_names(validate=False)
    except (DictionaryManifestError, FileNotFoundError):
        return frozenset()
    return frozenset(names)


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

    for section in ("activity", "assay", "cellline", "tissue", "testitem", "document", "target"):
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


def load_config(
    path: str | Path | None = None,
    cli_overrides: dict[str, Any] | None = None,
    *,
    base_path: Path | str | None = None,
    strict: bool = True,
    include_metadata: bool = False,
    cli_sources: TypingMapping[tuple[str, ...], str] | None = None,
) -> Config | tuple[Config, ConfigMetadata]:
    """Load configuration from *path* applying overrides."""

    try:
        data, resolved_path = load_yaml_config(path)
    except ConfigLoaderError as exc:
        raise ConfigError(str(exc)) from exc

    cli_path_map: dict[str, tuple[str, ...]] = {}
    source_map: dict[tuple[str, ...], ConfigSource] = {}
    detail = _format_source_detail(resolved_path)
    for path_tuple, _ in _iter_leaf_items(data):
        source_map[path_tuple] = ConfigSource("config", detail)

    local_path = resolved_path.with_name(
        f"{resolved_path.stem}.local{resolved_path.suffix}"
    )
    if local_path != resolved_path and local_path.exists():
        try:
            local_data, _ = load_yaml_config(local_path)
        except ConfigLoaderError as exc:
            raise ConfigError(str(exc)) from exc
        _merge_mapping(data, local_data)
        local_detail = _format_source_detail(local_path)
        for path_tuple, _ in _iter_leaf_items(local_data):
            source_map[path_tuple] = ConfigSource("config", local_detail)

    env_overrides, env_warnings = _apply_env_overrides(data)
    for path_tuple, env_key in env_overrides.items():
        source_map[path_tuple] = ConfigSource("env", env_key)
    _upgrade_legacy_config(data)

    if cli_overrides:
        for key, val in cli_overrides.items():
            parts = key.split(".")
            _set_by_path(data, parts, val)
            path_tuple = tuple(parts)
            detail = cli_sources.get(path_tuple) if cli_sources else None
            source_map[path_tuple] = ConfigSource("cli", detail)
            if detail:
                cli_path_map[detail] = path_tuple

    placeholder_base = _resolve_placeholder_base_path(base_path)
    data = _expand_config_placeholders(data, base_path=placeholder_base)

    base_dir = resolved_path.parent.resolve()
    _absolutise_config_paths(data, base_dir)

    unknown = _collect_unknown_keys(data, Config)
    ignored_unknown = [key for key in unknown if key in _OPTIONAL_UNKNOWN_KEYS]
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

    configure_rate_limiters(cfg)
    if not include_metadata:
        return cfg

    cfg_dict = _serialize_paths(cfg.to_dict())
    for path_tuple, _ in _iter_leaf_items(cfg_dict):
        source_map.setdefault(path_tuple, ConfigSource("default", None))
    snapshot = _build_snapshot(cfg_dict, source_map)
    metadata = ConfigMetadata(snapshot=snapshot, sources=source_map)
    if env_warnings:
        metadata.env_warnings.extend(env_warnings)
    if cli_path_map:
        metadata.cli_paths.update(cli_path_map)
    return cfg, metadata


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
