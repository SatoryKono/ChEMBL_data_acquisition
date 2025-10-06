"""Configuration loading utilities."""

from __future__ import annotations

import logging
import os
from pathlib import Path
from types import UnionType
from typing import Any, Callable, Iterator, Mapping, Sequence, Union, get_args, get_origin

import yaml
from pydantic import BaseModel, ValidationError
from pydantic_core import ErrorDetails

from ..common.log import logger
from ..utils.config import ConfigLoaderError, load_yaml_config
from .models import Config, ConfigMetadata, ConfigSource, _default_base_path

__all__ = [
    "ConfigError",
    "_absolutise_path_value",
    "_mask_secrets",
    "_serialize_paths",
    "build_alias_map",
    "load_config",
    "print_config",
]


class ConfigError(RuntimeError):
    """Raised when configuration loading fails."""


def _annotation_is_path(annotation: Any) -> bool:
    try:
        return issubclass(annotation, Path)
    except TypeError:
        return False


def _annotation_is_model(annotation: Any) -> bool:
    try:
        return issubclass(annotation, BaseModel)
    except TypeError:
        return False


def _collect_path_field_paths(
    model: type[BaseModel], prefix: tuple[str, ...] = ()
) -> set[tuple[str, ...]]:
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


_CONFIG_PATH_FIELDS: set[tuple[str, ...]] = _collect_path_field_paths(Config)


def _absolutise_path_value(value: Any, base_dir: Path) -> Any:
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


def _expand_config_placeholders(data: Any, *, base_path: Path) -> Any:
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
    if base_path is not None:
        candidate = Path(base_path).expanduser()
        if candidate.is_absolute():
            return candidate.resolve()
        return (Path.cwd() / candidate).resolve()
    return _default_base_path()


def _iter_leaf_items(
    data: Any, prefix: tuple[str, ...] = ()
) -> Iterator[tuple[tuple[str, ...], Any]]:
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


def _set_by_path(data: dict[str, Any], path: list[str], value: Any) -> None:
    cur: dict[str, Any] = data
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]  # type: ignore[assignment]
    cur[path[-1]] = value


def _coerce_integral_numbers(value: Any) -> Any:
    if isinstance(value, float):
        return int(value) if value.is_integer() else value
    if isinstance(value, list):
        return [_coerce_integral_numbers(item) for item in value]
    if isinstance(value, tuple):
        return tuple(_coerce_integral_numbers(item) for item in value)
    if isinstance(value, dict):
        return {key: _coerce_integral_numbers(val) for key, val in value.items()}
    return value


def _parse_env_value(env_key: str, raw_value: str) -> Any:
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


def _format_error_location(location: Sequence[Any]) -> str:
    parts: list[str] = []
    for part in location:
        if isinstance(part, str):
            parts.append(part)
        elif isinstance(part, int):
            parts.append(str(part))
    return ".".join(parts)


def _format_env_error_message(error: Mapping[str, Any]) -> str:
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
    if error_type == "bool_parsing":
        return "must be a boolean"
    if error_type == "float_parsing":
        return "must be a number"
    return error.get("msg", "is invalid")


def _format_env_error(env_key: str, error: Mapping[str, Any]) -> str:
    message = _format_env_error_message(error)
    loc = error.get("loc")
    location = _format_error_location(loc) if isinstance(loc, Sequence) else ""
    if location:
        return f"{env_key} ({location}) {message}".strip()
    return f"{env_key} {message}".strip()


def _normalize_env_errors(
    errors: Sequence[ErrorDetails], overrides: Mapping[tuple[str, ...], str]
) -> tuple[list[str], int]:
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


_OPTIONAL_UNKNOWN_KEYS: frozenset[str] = frozenset({
    "local.io.csv_fallback_separators",
})


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
                _collect_unknown_keys(val, submodel, prefix=f"{prefix}{key}."),
            )
    return unknown


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


def _merge_mapping(dest: dict[str, Any], src: dict[str, Any]) -> None:
    for key, value in src.items():
        if key in dest and isinstance(dest[key], dict) and isinstance(value, dict):
            _merge_mapping(dest[key], value)
        else:
            dest[key] = value


def _upgrade_legacy_config(data: dict[str, Any]) -> None:
    sources = data.setdefault("sources", {})
    chembl = sources.setdefault("chembl", {})
    local = data.setdefault("local", {})
    pipelines = chembl.setdefault("pipelines", {})
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

    for section in (
        "activity",
        "assay",
        "cellline",
        "tissue",
        "testitem",
        "document",
        "target",
    ):
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


def load_config(
    path: str | Path | None = None,
    cli_overrides: dict[str, Any] | None = None,
    *,
    base_path: Path | str | None = None,
    strict: bool = True,
    include_metadata: bool = False,
    cli_sources: Mapping[tuple[str, ...], str] | None = None,
) -> Config | tuple[Config, ConfigMetadata]:
    try:
        data, resolved_path = load_yaml_config(path)
    except ConfigLoaderError as exc:
        raise ConfigError(str(exc)) from exc

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

    data = _coerce_integral_numbers(data)

    if "$defs" in data:
        raise ConfigError(
            f"{resolved_path} appears to be a configuration schema; "
            "provide an application config file such as config/config.yaml."
        )

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

    from .runtime import configure_rate_limits

    configure_rate_limits(cfg.rate)

    if not include_metadata:
        return cfg

    cfg_dict = _serialize_paths(cfg.to_dict())
    for path_tuple, _ in _iter_leaf_items(cfg_dict):
        source_map.setdefault(path_tuple, ConfigSource("default", None))
    snapshot = _build_snapshot(cfg_dict, source_map)
    metadata = ConfigMetadata(snapshot=snapshot, sources=source_map)
    return cfg, metadata


def print_config(cfg: Config) -> None:
    data = _serialize_paths(cfg.to_dict())
    masked = _mask_secrets(data)
    print(yaml.safe_dump(masked, sort_keys=False))


def build_alias_map(
    model: type[BaseModel], prefix: str = "CHEMBL_DA"
) -> dict[str, list[str]]:
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
    "CHEMBL_DA_UNIPROT_BASE": ["sources", "uniprot", "base"],
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
