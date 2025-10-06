"""Helpers for loading and normalising configuration files."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, Mapping as TypingMapping

import yaml
from pydantic import BaseModel, ValidationError
from pydantic_core import ErrorDetails

from ..common.log import logger
from ..utils.config import ConfigLoaderError, load_yaml_config
from .models import (
    Config,
    ConfigMetadata,
    ConfigSource,
    ConfigError,
    _CONFIG_PATH_FIELDS,
    _ALIAS_MAP,
)
from .runtime import configure_rate_limiters

__all__ = [
    "ConfigMetadata",
    "ConfigSource",
    "ConfigError",
    "load_config",
    "ensure_dirs",
    "print_config",
    "_absolutise_path_value",
    "_serialize_paths",
    "_mask_secrets",
]


def _absolutise_path_value(value: Any, base_dir: Path) -> Any:
    """Return *value* converted to an absolute path relative to *base_dir*."""

    if value is None:
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


def _set_by_path(data: dict[str, Any], path: list[str], value: Any) -> None:
    cur: dict[str, Any] = data
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[path[-1]] = value


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
    """Normalize *raw_value* from environment variable *env_key*."""

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
    if error_type == "string_type":
        return "must be a string"
    if error_type == "bool_parsing":
        return "must be a boolean"
    if error_type == "finite_number":
        return "must be finite"
    return str(error.get("msg", "invalid value"))


def _format_env_error(env_key: str, error: Mapping[str, Any]) -> str:
    """Return a human readable error string for *env_key* based on *error*."""

    message = _format_env_error_message(error)
    loc = error.get("loc")
    location = _format_error_location(loc) if isinstance(loc, Sequence) else ""
    if location:
        return f"{env_key} ({location}) {message}".strip()
    return f"{env_key} {message}".strip()


def _normalize_env_errors(
    errors: Sequence[ErrorDetails], overrides: TypingMapping[tuple[str, ...], str]
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
    for path_tuple, _ in _iter_leaf_items(data):
        source_map[path_tuple] = ConfigSource("config", str(resolved_path))

    local_path = resolved_path.with_name(
        f"{resolved_path.stem}.local{resolved_path.suffix}"
    )
    if local_path != resolved_path and local_path.exists():
        try:
            local_data, _ = load_yaml_config(local_path)
        except ConfigLoaderError as exc:
            raise ConfigError(str(exc)) from exc
        _merge_mapping(data, local_data)
        for path_tuple, _ in _iter_leaf_items(local_data):
            source_map[path_tuple] = ConfigSource("config", str(local_path))

    env_overrides = _apply_env_overrides(data)
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
