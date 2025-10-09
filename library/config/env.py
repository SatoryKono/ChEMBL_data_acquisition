"""Environment utilities for configuration loading."""

from __future__ import annotations

import os
from collections.abc import Mapping as TypingMapping
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import yaml  # type: ignore[import-untyped]
from pydantic import BaseModel
from pydantic_core import ErrorDetails

from ..common.log import logger

if TYPE_CHECKING:  # pragma: no cover - imported only for typing
    from .models import Config

__all__ = [
    "_normalize_base_path",
    "_default_base_path",
    "_default_cache_home",
    "_resolve_placeholder_base_path",
    "_expand_config_placeholders",
    "_coerce_integral_numbers",
    "_parse_env_value",
    "_format_env_error_message",
    "_format_env_error",
    "_normalize_env_errors",
    "_apply_env_overrides",
    "_set_by_path",
]


def _config_model() -> type[Config]:
    from .models import Config

    return Config


def _alias_map() -> dict[str, tuple[str, ...]]:
    from .models import _ALIAS_MAP

    # Ensure all values are tuples to match the expected return type
    return {k: tuple(v) for k, v in _ALIAS_MAP.items()}

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


def _resolve_placeholder_base_path(base_path: Path | str | None) -> Path:
    """Return the base path used for ``$CHEMBL_DA_BASE_PATH`` placeholders."""

    if base_path is not None:
        return _normalize_base_path(base_path)
    return _default_base_path()


def _expand_config_placeholders(data: Any, *, base_path: Path) -> Any:
    """Expand configuration placeholders in ``data`` using *base_path*."""

    # Normalise to platform-specific path string
    base_str = str(Path(base_path))
    replacements = {"$CHEMBL_DA_BASE_PATH": base_str}

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
            had_placeholder = False
            for marker, target in replacements.items():
                if marker in replaced:
                    had_placeholder = True
                    replaced = replaced.replace(marker, target)
            replaced = os.path.expandvars(replaced)
            replaced = os.path.expanduser(replaced)
            # Only normalise path separators when placeholder was present
            if had_placeholder or replaced.startswith(base_str):
                try:
                    replaced = str(Path(replaced))
                except Exception:
                    return replaced
            return replaced
        return value

    return _expand(data)


def _set_by_path(data: dict[str, Any], path: Sequence[str], value: Any) -> None:
    cur: dict[str, Any] = data
    for key in path[:-1]:
        if key not in cur or not isinstance(cur[key], dict):
            cur[key] = {}
        cur = cur[key]
    cur[path[-1]] = value


def _is_valid_path(path: Sequence[str]) -> bool:
    model: type[BaseModel] | None = _config_model()
    for part in path:
        if model is None or part not in model.model_fields:
            return False
        field = model.model_fields[part].annotation
        model = (
            field if isinstance(field, type) and issubclass(field, BaseModel) else None
        )
    return True


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


def _parse_env_value(env_key: str, raw_value: str) -> Any:
    """Normalise *raw_value* from environment variable *env_key*."""

    if raw_value == "":
        return ""
    if raw_value and raw_value.strip() == "":
        return raw_value
    try:
        value = yaml.safe_load(raw_value)
    except yaml.YAMLError as exc:  # pragma: no cover - defensive logging
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


def _format_env_error_message(error: TypingMapping[str, Any]) -> str:
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


def _format_env_error(env_key: str, error: TypingMapping[str, Any]) -> str:
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
        handled += 1
        messages.append(_format_env_error(env_key, error))
    return messages, handled


def _apply_env_overrides(data: dict[str, Any]) -> dict[tuple[str, ...], str]:
    prefix = "CHEMBL_DA"
    overrides: dict[tuple[str, ...], str] = {}
    alias_map = _alias_map()
    for env_key, env_val in os.environ.items():
        key = env_key.upper()
        if key in alias_map:
            path: tuple[str, ...]
            if key in alias_map:
                path = alias_map[key]
            elif key.startswith(prefix + "__"):
                path = tuple(key[len(prefix) + 2 :].split("__"))
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

