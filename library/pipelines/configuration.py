"""Loader for declarative pipeline metadata stored under :mod:`config.pipeline`."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from typing import Any

import yaml

from config.paths import PIPELINE_DIR
from library.pipelines.common.metadata import get_pipeline_version

__all__ = [
    "PipelineConfig",
    "PipelineIO",
    "PipelineParameters",
    "PipelineStepConfig",
    "PipelineVersionInfo",
    "list_pipeline_configs",
    "load_pipeline_config",
]


@dataclass(frozen=True)
class PipelineVersionInfo:
    """Metadata describing the pipeline version column used in exports."""

    column: str
    value: str
    description: str | None = None


@dataclass(frozen=True)
class PipelineStepConfig:
    """Declarative description of a pipeline step."""

    name: str
    description: str | None
    enabled: bool
    depends_on: tuple[str, ...]
    produces: tuple[str, ...]
    applies_to: tuple[str, ...]


@dataclass(frozen=True)
class PipelineParameters:
    """Structured mapping of CLI parameters to configuration paths."""

    default: Mapping[str, str]
    shared: Mapping[str, str]
    modes: Mapping[str, Mapping[str, str]]
    commands: Mapping[str, Mapping[str, str]]

    def mapping_for(
        self, *, mode: str | None = None, command: str | None = None
    ) -> dict[str, str]:
        """Return a merged mapping using ``mode``/``command`` specific overrides."""

        mapping: dict[str, str] = dict(self.default)
        if self.shared:
            mapping.update(self.shared)
        if mode is not None and mode in self.modes:
            mapping.update(self.modes[mode])
        if command is not None and command in self.commands:
            mapping.update(self.commands[command])
        return mapping


@dataclass(frozen=True)
class PipelineIO:
    """Inputs or outputs declared for a pipeline."""

    primary: str | None
    dependencies: tuple[str, ...]
    variants: Mapping[str, Any]


@dataclass(frozen=True)
class PipelineConfig:
    """Complete declarative definition of a pipeline."""

    name: str
    summary: str | None
    schemas: Mapping[str, str]
    pipeline_version: PipelineVersionInfo
    inputs: PipelineIO
    outputs: PipelineIO
    parameters: PipelineParameters
    steps: tuple[PipelineStepConfig, ...]


def _coerce_str(value: Any, *, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise TypeError(f"pipeline config field '{field}' must be a non-empty string")
    return value.strip()


def _coerce_optional_str(value: Any, *, field: str) -> str | None:
    if value is None:
        return None
    return _coerce_str(value, field=field)


def _coerce_str_sequence(value: Any, *, field: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (_coerce_str(value, field=field),)
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        items: list[str] = []
        for idx, item in enumerate(value):
            items.append(_coerce_str(item, field=f"{field}[{idx}]"))
        return tuple(items)
    raise TypeError(
        f"pipeline config field '{field}' must be a string or sequence of strings"
    )


def _coerce_mapping_str_str(value: Any, *, field: str) -> dict[str, str]:
    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise TypeError(f"pipeline config field '{field}' must be a mapping")
    result: dict[str, str] = {}
    for key, val in value.items():
        result[_coerce_str(key, field=f"{field} key")] = _coerce_str(
            val, field=f"{field}[{key}]"
        )
    return result


def _parse_pipeline_version(
    data: Mapping[str, Any], *, name: str
) -> PipelineVersionInfo:
    raw = data.get("pipeline_version")
    if raw is None:
        raise KeyError(f"pipeline '{name}' is missing 'pipeline_version'")
    column: str
    description: str | None = None
    value: str | None = None
    if isinstance(raw, str):
        column = _coerce_str(raw, field="pipeline_version")
    elif isinstance(raw, Mapping):
        column = _coerce_str(raw.get("column"), field="pipeline_version.column")
        description = _coerce_optional_str(
            raw.get("description"), field="pipeline_version.description"
        )
        raw_value = raw.get("value")
        if raw_value is not None:
            value = _coerce_str(raw_value, field="pipeline_version.value")
    else:
        raise TypeError("pipeline_version must be a string or mapping")
    if value is None or value.lower() == "auto":
        value = get_pipeline_version()
    return PipelineVersionInfo(column=column, value=value, description=description)


def _parse_steps(data: Mapping[str, Any]) -> tuple[PipelineStepConfig, ...]:
    raw_steps = data.get("steps")
    if raw_steps is None:
        return ()
    if not isinstance(raw_steps, Sequence) or isinstance(
        raw_steps, (str, bytes, bytearray)
    ):
        raise TypeError("pipeline 'steps' must be a sequence of mappings")
    steps: list[PipelineStepConfig] = []
    for index, entry in enumerate(raw_steps):
        if not isinstance(entry, Mapping):
            raise TypeError(f"pipeline step at index {index} must be a mapping")
        name = _coerce_str(entry.get("name"), field=f"steps[{index}].name")
        description = _coerce_optional_str(
            entry.get("description"), field=f"steps[{index}].description"
        )
        enabled_raw = entry.get("enabled", True)
        if enabled_raw in (None, "", "auto"):
            enabled = True
        elif isinstance(enabled_raw, bool):
            enabled = enabled_raw
        else:
            enabled = str(enabled_raw).lower() not in {"false", "0", "no", "off"}
        depends_on = _coerce_str_sequence(
            entry.get("depends_on"), field=f"steps[{index}].depends_on"
        )
        produces = _coerce_str_sequence(
            entry.get("produces"), field=f"steps[{index}].produces"
        )
        applies_to = _coerce_str_sequence(
            entry.get("applies_to"), field=f"steps[{index}].applies_to"
        )
        steps.append(
            PipelineStepConfig(
                name=name,
                description=description,
                enabled=enabled,
                depends_on=depends_on,
                produces=produces,
                applies_to=applies_to,
            )
        )
    return tuple(steps)


def _parse_parameters(data: Mapping[str, Any]) -> PipelineParameters:
    raw_params = data.get("parameters")
    if raw_params is None:
        return PipelineParameters(default={}, shared={}, modes={}, commands={})
    if not isinstance(raw_params, Mapping):
        raise TypeError("pipeline 'parameters' must be a mapping")
    default: dict[str, str] = {}
    shared: dict[str, str] = {}
    modes: dict[str, dict[str, str]] = {}
    commands: dict[str, dict[str, str]] = {}

    def _consume_plain_entries(
        entries: Mapping[str, Any], *, dest: dict[str, str], prefix: str
    ) -> None:
        for key, value in entries.items():
            if key in {"shared", "modes", "commands", "default"}:
                continue
            if isinstance(value, Mapping):
                raise TypeError(
                    "nested mappings within 'parameters' must be under the 'shared', 'modes', "
                    "'commands', or 'default' keys"
                )
            dest[_coerce_str(key, field=f"parameters.{prefix}{key}")] = _coerce_str(
                value, field=f"parameters.{prefix}{key}"
            )

    if "default" in raw_params:
        default.update(
            _coerce_mapping_str_str(
                raw_params.get("default"), field="parameters.default"
            )
        )
    _consume_plain_entries(raw_params, dest=default, prefix="")
    shared.update(
        _coerce_mapping_str_str(raw_params.get("shared"), field="parameters.shared")
    )

    raw_modes = raw_params.get("modes")
    if raw_modes is not None:
        if not isinstance(raw_modes, Mapping):
            raise TypeError("parameters.modes must be a mapping")
        for mode, mapping in raw_modes.items():
            modes[_coerce_str(mode, field="parameters.modes key")] = (
                _coerce_mapping_str_str(mapping, field=f"parameters.modes.{mode}")
            )

    raw_commands = raw_params.get("commands")
    if raw_commands is not None:
        if not isinstance(raw_commands, Mapping):
            raise TypeError("parameters.commands must be a mapping")
        for command, mapping in raw_commands.items():
            commands[_coerce_str(command, field="parameters.commands key")] = (
                _coerce_mapping_str_str(mapping, field=f"parameters.commands.{command}")
            )

    return PipelineParameters(
        default=default, shared=shared, modes=modes, commands=commands
    )


def _parse_io_section(data: Mapping[str, Any], *, field: str) -> PipelineIO:
    raw_io = data.get(field)
    if raw_io is None:
        return PipelineIO(primary=None, dependencies=(), variants={})
    if not isinstance(raw_io, Mapping):
        raise TypeError(f"pipeline '{field}' must be a mapping")
    primary = raw_io.get("primary")
    dependencies = raw_io.get("dependencies")
    variants = {k: v for k, v in raw_io.items() if k not in {"primary", "dependencies"}}
    primary_str = _coerce_optional_str(primary, field=f"{field}.primary")
    deps = _coerce_str_sequence(dependencies, field=f"{field}.dependencies")
    return PipelineIO(primary=primary_str, dependencies=deps, variants=variants)


def _parse_schemas(data: Mapping[str, Any], *, name: str) -> dict[str, str]:
    raw = data.get("schemas")
    if raw is None:
        raise KeyError(f"pipeline '{name}' is missing the 'schemas' section")
    return _coerce_mapping_str_str(raw, field="schemas")


def _ensure_directory(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_dir():
        raise FileNotFoundError(
            f"pipeline configuration directory not found: {resolved}"
        )
    return resolved


@cache
def list_pipeline_configs(directory: Path | None = None) -> tuple[str, ...]:
    """Return the available pipeline configuration names."""

    base = _ensure_directory(Path(directory) if directory else PIPELINE_DIR)
    names = sorted(path.stem for path in base.glob("*.yaml"))
    return tuple(names)


@cache
def load_pipeline_config(name: str, *, directory: Path | None = None) -> PipelineConfig:
    """Return the parsed declarative configuration for ``name``."""

    if not name:
        raise ValueError("pipeline name must be provided")
    base = _ensure_directory(Path(directory) if directory else PIPELINE_DIR)
    path = base / f"{name}.yaml"
    if not path.exists():
        available = ", ".join(list_pipeline_configs(directory=base))
        raise FileNotFoundError(
            f"pipeline configuration not found: {path} (available: {available})"
        )
    with path.open("r", encoding="utf-8") as handle:
        raw_data = yaml.safe_load(handle) or {}
    if not isinstance(raw_data, Mapping):
        raise TypeError(f"pipeline configuration '{path}' must contain a mapping")
    cfg_name = _coerce_str(raw_data.get("name", name), field="name")
    summary = _coerce_optional_str(raw_data.get("summary"), field="summary")
    schemas = _parse_schemas(raw_data, name=cfg_name)
    pipeline_version = _parse_pipeline_version(raw_data, name=cfg_name)
    parameters = _parse_parameters(raw_data)
    steps = _parse_steps(raw_data)
    inputs = _parse_io_section(raw_data, field="inputs")
    outputs = _parse_io_section(raw_data, field="outputs")
    return PipelineConfig(
        name=cfg_name,
        summary=summary,
        schemas=schemas,
        pipeline_version=pipeline_version,
        inputs=inputs,
        outputs=outputs,
        parameters=parameters,
        steps=steps,
    )
