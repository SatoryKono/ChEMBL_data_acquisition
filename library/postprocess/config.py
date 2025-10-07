"""Load and normalise declarative post-processing pipeline definitions."""
from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Iterator, Mapping, Sequence

import yaml

from config.paths import PIPELINE_DIR as _PIPELINE_DIR
from library.config.env import _expand_config_placeholders, _resolve_placeholder_base_path

__all__ = [
    "PIPELINE_CONFIG_DIR",
    "PipelineConfig",
    "PipelineConfigError",
    "PipelineStep",
    "load_pipeline_config",
]

PIPELINE_CONFIG_DIR = _PIPELINE_DIR


class PipelineConfigError(RuntimeError):
    """Raised when a pipeline configuration file cannot be parsed."""


@dataclass(frozen=True)
class PipelineStep:
    """Enabled step definition declared in a pipeline configuration file."""

    name: str
    callable_path: str
    params: Mapping[str, Any]

    def resolve(self) -> Callable[..., Any]:
        """Return the callable referenced by :attr:`callable_path`."""

        module_name, _, attribute = self.callable_path.partition(":")
        if not module_name or not attribute:
            raise PipelineConfigError(
                f"invalid callable reference '{self.callable_path}' for step '{self.name}'"
            )
        try:
            module = import_module(module_name)
        except ImportError as exc:  # pragma: no cover - defensive
            raise PipelineConfigError(
                f"module '{module_name}' for step '{self.name}' cannot be imported"
            ) from exc
        try:
            return getattr(module, attribute)
        except AttributeError as exc:  # pragma: no cover - defensive
            raise PipelineConfigError(
                f"callable '{self.callable_path}' for step '{self.name}' is not importable"
            ) from exc


@dataclass(frozen=True)
class PipelineConfig:
    """Top-level representation of a pipeline configuration file."""

    pipeline_version: str
    steps: tuple[PipelineStep, ...]
    path: Path

    def iter_callables(self) -> Iterator[Callable[..., Any]]:
        """Yield callables for enabled steps in declaration order."""

        for step in self.steps:
            yield step.resolve()


def load_pipeline_config(
    domain: str,
    *,
    path: str | Path | None = None,
    base_path: str | Path | None = None,
) -> PipelineConfig:
    """Return post-processing pipeline configuration for ``domain``."""

    config_path = _resolve_config_path(domain, path)
    data = _load_yaml(config_path)
    expanded = _expand_config_placeholders(
        data,
        base_path=_resolve_placeholder_base_path(base_path),
    )
    pipeline_version = expanded.get("pipeline_version")
    if not isinstance(pipeline_version, str) or not pipeline_version.strip():
        raise PipelineConfigError(
            f"pipeline configuration '{config_path}' is missing a non-empty 'pipeline_version'"
        )

    raw_steps = expanded.get("steps", [])
    if not isinstance(raw_steps, Sequence):
        raise PipelineConfigError(
            f"pipeline configuration '{config_path}' must define 'steps' as a sequence"
        )

    steps: list[PipelineStep] = []
    for entry in raw_steps:
        if not isinstance(entry, Mapping):
            raise PipelineConfigError(
                f"pipeline configuration '{config_path}' contains a non-mapping step entry"
            )
        enabled = entry.get("enabled", True)
        if not enabled:
            continue
        name = entry.get("name")
        callable_path = entry.get("callable")
        params = entry.get("params", {})
        if not isinstance(name, str) or not name.strip():
            raise PipelineConfigError(
                f"pipeline configuration '{config_path}' has an enabled step without a string 'name'"
            )
        if not isinstance(callable_path, str) or not callable_path.strip():
            raise PipelineConfigError(
                f"pipeline configuration '{config_path}' step '{name}' is missing a 'callable'"
            )
        if params is None:
            params = {}
        if not isinstance(params, Mapping):
            raise PipelineConfigError(
                f"pipeline configuration '{config_path}' step '{name}' has non-mapping 'params'"
            )
        steps.append(
            PipelineStep(
                name=name,
                callable_path=callable_path,
                params=_freeze_params(params),
            )
        )

    return PipelineConfig(
        pipeline_version=pipeline_version,
        steps=tuple(steps),
        path=config_path,
    )


def _resolve_config_path(domain: str, path: str | Path | None) -> Path:
    if path is not None:
        return Path(path).expanduser().resolve()
    candidate = PIPELINE_CONFIG_DIR / f"{domain}.yaml"
    return candidate.resolve()


def _load_yaml(path: Path) -> Mapping[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    except FileNotFoundError as exc:
        raise PipelineConfigError(f"pipeline configuration file not found: {path}") from exc
    except yaml.YAMLError as exc:
        raise PipelineConfigError(f"failed to parse YAML at {path}: {exc}") from exc
    if not isinstance(data, Mapping):
        raise PipelineConfigError(
            f"pipeline configuration '{path}' must contain a top-level mapping"
        )
    return data


def _freeze_params(params: Mapping[str, Any]) -> Mapping[str, Any]:
    def _convert(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {str(key): _convert(val) for key, val in value.items()}
        if isinstance(value, list):
            return [_convert(item) for item in value]
        if isinstance(value, tuple):
            return tuple(_convert(item) for item in value)
        if isinstance(value, set):
            return sorted(_convert(item) for item in value)
        return value

    materialised = {str(key): _convert(val) for key, val in params.items()}
    return MappingProxyType(materialised)
