"""Load post-processing pipeline configurations from YAML files."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping

import os
import re

import yaml

from config.paths import PIPELINE_CONFIG_DIR

from .common.import_utils import resolve_dotted_path
from .common.types import StepDefinition

__all__ = ["PipelineConfig", "PipelineStep", "load_pipeline_config"]


_ENV_PATTERN = re.compile(r"\$\{([^:{}\s]+)(:-([^{}]*))?\}")


@dataclass(frozen=True)
class PipelineStep:
    """Declarative description of a configured pipeline step."""

    name: str
    dotted_path: str
    enabled: bool
    params: Mapping[str, Any]
    description: str | None = None


@dataclass(frozen=True)
class PipelineConfig:
    """Fully resolved pipeline configuration."""

    name: str
    pipeline_version: str
    steps: tuple[StepDefinition, ...]
    raw_steps: tuple[PipelineStep, ...]


def _interpolate_environment(raw: str, env: Mapping[str, str]) -> str:
    """Expand ``${VAR}`` or ``${VAR:-default}`` placeholders in *raw*."""

    def _replace(match: re.Match[str]) -> str:
        key = match.group(1)
        default = match.group(3)
        value = env.get(key)
        if value is None:
            return default or ""
        return value

    expanded = _ENV_PATTERN.sub(_replace, raw)
    expanded = os.path.expandvars(expanded)
    expanded = os.path.expanduser(expanded)
    return expanded


def _load_yaml(path: Path, env: Mapping[str, str]) -> Mapping[str, Any]:
    """Return parsed YAML from *path* with environment interpolation."""

    raw = path.read_text(encoding="utf-8")
    interpolated = _interpolate_environment(raw, env)
    data = yaml.safe_load(interpolated) or {}
    if not isinstance(data, Mapping):  # pragma: no cover - defensive guard
        raise TypeError(f"pipeline config at {path} must be a mapping")
    return data


def _coerce_mapping(value: Any, *, context: str) -> Mapping[str, Any]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError(f"{context} must be a mapping")
    return value


def _build_step(config: PipelineStep) -> StepDefinition:
    func = resolve_dotted_path(config.dotted_path)

    if config.params:

        def _wrapped(df, _func=func, _params=dict(config.params)):
            return _func(df, **_params)

        step_func = _wrapped
    else:

        def _wrapped(df, _func=func):
            return _func(df)

        step_func = _wrapped

    return StepDefinition(config.name, step_func, config.description)


def _parse_steps(data: Iterable[Any], *, name: str) -> tuple[PipelineStep, ...]:
    steps: list[PipelineStep] = []
    for index, entry in enumerate(data, start=1):
        if not isinstance(entry, Mapping):
            raise TypeError(f"step #{index} in pipeline '{name}' must be a mapping")
        raw_name = entry.get("name")
        if raw_name is None:
            raise ValueError(f"step #{index} in pipeline '{name}' missing 'name'")
        step_name = str(raw_name)
        if not step_name.strip():
            raise ValueError(f"step #{index} in pipeline '{name}' missing 'name'")
        dotted = entry.get("callable")
        if not isinstance(dotted, str) or not dotted:
            raise ValueError(
                f"step '{step_name}' in pipeline '{name}' missing callable"
            )
        enabled = bool(entry.get("enabled", True))
        params = _coerce_mapping(entry.get("params"), context=f"step '{step_name}' params")
        description_value = entry.get("description")
        description = str(description_value) if description_value is not None else None
        steps.append(
            PipelineStep(
                name=step_name,
                dotted_path=dotted,
                enabled=enabled,
                params=params,
                description=description,
            )
        )
    return tuple(steps)


def _load_pipeline(path: Path, *, env: Mapping[str, str]) -> PipelineConfig:
    data = _load_yaml(path, env)
    pipeline_version = str(data.get("pipeline_version", "")) or "dev"
    raw_steps_data = data.get("steps") or []
    if not isinstance(raw_steps_data, list):
        raise TypeError(f"pipeline '{path.stem}' steps must be a list")
    parsed_steps = _parse_steps(raw_steps_data, name=path.stem)
    active_steps = [step for step in parsed_steps if step.enabled]
    definitions = tuple(_build_step(step) for step in active_steps)
    return PipelineConfig(
        name=path.stem,
        pipeline_version=pipeline_version,
        steps=definitions,
        raw_steps=parsed_steps,
    )


@lru_cache(maxsize=None)
def _load_pipeline_cached(name: str, base_dir: str | None) -> PipelineConfig:
    base = Path(base_dir) if base_dir is not None else PIPELINE_CONFIG_DIR
    path = (base / f"{name}.yaml").resolve()
    if not path.exists():  # pragma: no cover - defensive guard
        raise FileNotFoundError(f"pipeline configuration not found: {path}")
    return _load_pipeline(path, env=dict(os.environ))


def load_pipeline_config(
    name: str,
    *,
    base_dir: Path | None = None,
    env: Mapping[str, str] | None = None,
) -> PipelineConfig:
    """Return the resolved pipeline configuration for *name*."""

    if env is None:
        base_str = None if base_dir is None else str(base_dir)
        return _load_pipeline_cached(name, base_str)

    base = base_dir or PIPELINE_CONFIG_DIR
    path = (base / f"{name}.yaml").resolve()
    if not path.exists():  # pragma: no cover - defensive guard
        raise FileNotFoundError(f"pipeline configuration not found: {path}")
    return _load_pipeline(path, env=dict(env))
