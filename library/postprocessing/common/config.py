"""Helpers for loading declarative post-processing pipeline configuration."""
from __future__ import annotations

import inspect
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml

from config.paths import PIPELINE_DIR

from .import_utils import import_by_path
from .types import StepDefinition

__all__ = [
    "ConfiguredStep",
    "PipelineConfig",
    "PipelineConfigError",
    "load_pipeline_config",
    "normalize_pipeline_version",
]


_ENV_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)(?::(-|:-)([^}]*))?\}")


class PipelineConfigError(RuntimeError):
    """Raised when a pipeline configuration file cannot be loaded."""


@dataclass(frozen=True)
class ConfiguredStep:
    """A post-processing step loaded from the declarative configuration."""

    name: str
    callable_path: str
    definition: StepDefinition
    params: Mapping[str, Any]
    description: str | None = None


@dataclass(frozen=True)
class PipelineConfig:
    """Container holding the resolved pipeline configuration."""

    name: str
    path: Path
    pipeline_version: str | None
    params: Mapping[str, Any]
    steps: tuple[ConfiguredStep, ...]

    def step_definitions(self) -> tuple[StepDefinition, ...]:
        """Return the resolved step definitions preserving configuration order."""

        return tuple(step.definition for step in self.steps)


def load_pipeline_config(
    name: str,
    path: Path | None = None,
    *,
    env: Mapping[str, str] | None = None,
) -> PipelineConfig:
    """Load the pipeline configuration identified by ``name``.

    Parameters
    ----------
    name:
        Domain identifier, for example ``"activities"`` or ``"targets"``.
    path:
        Optional path overriding the default ``config/pipeline/<name>.yaml`` location.
    env:
        Mapping used to resolve environment variable placeholders. Defaults to
        :data:`os.environ` when omitted.
    """

    config_path = _resolve_config_path(name, path)
    env_mapping = dict(os.environ if env is None else env)

    try:
        raw_data = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    except FileNotFoundError as exc:  # pragma: no cover - defensive guard
        raise PipelineConfigError(
            f"pipeline configuration not found: {config_path}"
        ) from exc
    except yaml.YAMLError as exc:  # pragma: no cover - defensive guard
        raise PipelineConfigError(
            f"failed to parse pipeline configuration at {config_path}: {exc}"
        ) from exc

    if not isinstance(raw_data, dict):
        raise PipelineConfigError("pipeline configuration must be a mapping at the top level")

    expanded = _expand_environment(raw_data, env_mapping)

    pipeline_version = expanded.get("pipeline_version")
    if pipeline_version is not None and not isinstance(pipeline_version, str):
        raise PipelineConfigError("'pipeline_version' must be a string when provided")

    params = expanded.get("params") or {}
    if not isinstance(params, dict):
        raise PipelineConfigError("'params' must be a mapping when provided")

    enabled_steps = expanded.get("enabled_steps")
    if not isinstance(enabled_steps, Sequence) or isinstance(enabled_steps, (str, bytes)):
        raise PipelineConfigError("'enabled_steps' must be a sequence of mappings")

    steps = tuple(_load_step(entry, index) for index, entry in enumerate(enabled_steps))

    return PipelineConfig(
        name=name,
        path=config_path,
        pipeline_version=pipeline_version,
        params=params,
        steps=steps,
    )


def normalize_pipeline_version(value: str | None) -> str | None:
    """Return ``value`` normalised for runtime consumption.

    ``None``, empty strings and the sentinel values ``"auto"`` / ``"default"``
    become :data:`None`, signalling that callers should fall back to the
    library-reported pipeline version.
    """

    if value is None:
        return None
    candidate = value.strip()
    if not candidate:
        return None
    if candidate.lower() in {"auto", "default"}:
        return None
    return candidate


def _resolve_config_path(name: str, override: Path | None) -> Path:
    if override is not None:
        return Path(override).expanduser().resolve()
    return (PIPELINE_DIR / f"{name}.yaml").resolve()


def _load_step(entry: Any, index: int) -> ConfiguredStep:
    if not isinstance(entry, dict):
        raise PipelineConfigError(f"step #{index} must be a mapping")

    name = entry.get("name")
    callable_path = entry.get("callable")
    description = entry.get("description")
    params = entry.get("params") or {}

    if not isinstance(name, str) or not name:
        raise PipelineConfigError(f"step #{index} missing a valid 'name' string")
    if not isinstance(callable_path, str) or not callable_path:
        raise PipelineConfigError(f"step '{name}' missing a 'callable' reference")
    if description is not None and not isinstance(description, str):
        raise PipelineConfigError(f"step '{name}' description must be a string when provided")
    if not isinstance(params, dict):
        raise PipelineConfigError(f"step '{name}' parameters must be expressed as a mapping")

    func = import_by_path(callable_path)
    if not callable(func):  # pragma: no cover - defensive guard
        raise PipelineConfigError(
            f"resolved object for step '{name}' is not callable: {callable_path}"
        )

    _validate_step_parameters(func, params, step_name=name)

    definition = StepDefinition(
        name=name,
        func=func,
        description=description,
        params=params,
    )
    return ConfiguredStep(
        name=name,
        callable_path=callable_path,
        definition=definition,
        params=params,
        description=description,
    )


def _validate_step_parameters(func: Any, params: Mapping[str, Any], *, step_name: str) -> None:
    """Ensure that every configured parameter is supported by ``func``."""

    if not params:
        return

    try:
        signature = inspect.signature(func)
    except (TypeError, ValueError):  # pragma: no cover - defensive guard
        return

    if any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    ):
        return

    unsupported = sorted(
        key for key in params.keys() if key not in signature.parameters
    )
    if unsupported:
        formatted = ", ".join(unsupported)
        raise PipelineConfigError(
            f"step '{step_name}' defines unsupported parameters: {formatted}"
        )


def _expand_environment(value: Any, env: Mapping[str, str]) -> Any:
    if isinstance(value, dict):
        return {key: _expand_environment(item, env) for key, item in value.items()}
    if isinstance(value, list):
        return [_expand_environment(item, env) for item in value]
    if isinstance(value, tuple):  # pragma: no cover - handled implicitly via lists
        return tuple(_expand_environment(item, env) for item in value)
    if isinstance(value, str):
        return _expand_string(value, env)
    return value


def _expand_string(template: str, env: Mapping[str, str]) -> str:
    def _replacement(match: re.Match[str]) -> str:
        name = match.group(1)
        operator = match.group(2)
        default = match.group(3) or ""
        current = env.get(name)
        if operator == ":-":
            return current if current not in (None, "") else default
        if operator == "-":
            return current if current is not None else default
        return current or ""

    expanded = _ENV_VAR_PATTERN.sub(_replacement, template)
    expanded = os.path.expandvars(expanded)
    expanded = os.path.expanduser(expanded)
    return expanded
