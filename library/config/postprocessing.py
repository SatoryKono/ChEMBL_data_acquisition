"""Loader and validation helpers for post-processing pipeline configuration."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Sequence

import json
import sys

import yaml
from jsonschema import Draft202012Validator, FormatChecker
from jsonschema.exceptions import ValidationError as JsonSchemaValidationError

from config.paths import CONFIG_SCHEMA_PATH, POSTPROCESSING_CONFIG_DIR

__all__ = [
    "POSTPROCESSING_CONFIG_DIR",
    "PostprocessingConfigError",
    "PostprocessingPipeline",
    "PostprocessingStep",
    "diff_postprocessing_pipeline",
    "list_postprocessing_tables",
    "load_postprocessing_pipeline",
    "main",
    "validate_postprocessing_config",
]


SCHEMA_DEFINITION_REF = "#/$defs/PostprocessPipelineCfg"
_STEP_DEFAULT_KEYS = {
    "enabled",
    "continue_on_error",
    "kwargs",
    "gates",
    "description",
}


class PostprocessingConfigError(RuntimeError):
    """Raised when post-processing configuration parsing fails."""


@dataclass(frozen=True, slots=True)
class PostprocessingStep:
    """Description of a post-processing step declared in the pipeline config."""

    identifier: str
    callable: str
    enabled: bool = True
    continue_on_error: bool = False
    description: str | None = None
    gates: tuple[str, ...] = field(default_factory=tuple)
    kwargs: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:  # pragma: no cover - dataclass guard
        object.__setattr__(self, "gates", tuple(self.gates))
        object.__setattr__(self, "kwargs", MappingProxyType(dict(self.kwargs)))

    def as_dict(self) -> dict[str, Any]:
        """Return a serialisable representation of the step."""

        payload: dict[str, Any] = {
            "id": self.identifier,
            "callable": self.callable,
            "enabled": self.enabled,
            "continue_on_error": self.continue_on_error,
        }
        if self.description is not None:
            payload["description"] = self.description
        if self.gates:
            payload["gates"] = list(self.gates)
        if self.kwargs:
            payload["kwargs"] = {key: self.kwargs[key] for key in sorted(self.kwargs)}
        return payload


@dataclass(frozen=True, slots=True)
class PostprocessingPipeline:
    """Merged configuration describing the post-processing pipeline."""

    version: int
    flags: Mapping[str, bool]
    steps: tuple[PostprocessingStep, ...]
    step_defaults: Mapping[str, Any]

    def __post_init__(self) -> None:  # pragma: no cover - dataclass guard
        object.__setattr__(self, "flags", MappingProxyType(dict(self.flags)))
        object.__setattr__(self, "steps", tuple(self.steps))
        object.__setattr__(self, "step_defaults", MappingProxyType(dict(self.step_defaults)))

    def as_dict(self) -> dict[str, Any]:
        """Return a JSON/YAML compatible representation of the pipeline."""

        payload: dict[str, Any] = {
            "version": self.version,
            "flags": dict(self.flags),
            "steps": [step.as_dict() for step in self.steps],
        }
        if self.step_defaults:
            defaults: dict[str, Any] = {}
            for key, value in self.step_defaults.items():
                if key == "kwargs" and isinstance(value, Mapping):
                    defaults[key] = {k: value[k] for k in sorted(value)}
                else:
                    defaults[key] = value
            payload["defaults"] = {"step": defaults}
        return payload


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    """Return YAML content from *path* ensuring the top-level is a mapping."""

    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:  # pragma: no cover - defensive I/O guard
        raise PostprocessingConfigError(f"unable to read config: {path!s}") from exc

    try:
        data = yaml.safe_load(text) or {}
    except yaml.YAMLError as exc:
        raise PostprocessingConfigError(
            f"failed to parse YAML configuration at {path!s}: {exc}"
        ) from exc

    if not isinstance(data, dict):
        raise PostprocessingConfigError(
            f"top-level structure in {path!s} must be a mapping"
        )
    return data


@lru_cache(maxsize=None)
def _load_schema(schema_path: Path) -> Draft202012Validator:
    """Return the JSON schema validator for the post-processing config."""

    try:
        schema_text = schema_path.read_text(encoding="utf-8")
    except OSError as exc:  # pragma: no cover - defensive I/O guard
        raise PostprocessingConfigError(
            f"unable to read JSON schema: {schema_path!s}"
        ) from exc

    try:
        schema_data = json.loads(schema_text)
    except json.JSONDecodeError as exc:
        raise PostprocessingConfigError(
            f"invalid JSON schema at {schema_path!s}: {exc}"
        ) from exc

    scope: dict[str, Any] = {"$ref": SCHEMA_DEFINITION_REF}
    defs = schema_data.get("$defs")
    if isinstance(defs, dict):
        scope["$defs"] = defs
    return Draft202012Validator(scope, format_checker=FormatChecker())


def _validate_schema(instance: Mapping[str, Any], schema_path: Path) -> None:
    """Validate *instance* against the post-processing schema."""

    validator = _load_schema(schema_path)
    try:
        validator.validate(instance)
    except JsonSchemaValidationError as exc:  # pragma: no cover - schema guard
        raise PostprocessingConfigError(str(exc)) from exc


def _merge_step_defaults(
    base: Mapping[str, Any], override: Mapping[str, Any]
) -> dict[str, Any]:
    """Return combined defaults for a step."""

    merged: dict[str, Any] = {}
    for source in (base, override):
        defaults_section = source.get("defaults")
        if not isinstance(defaults_section, Mapping):
            continue
        step_section = defaults_section.get("step")
        if not isinstance(step_section, Mapping):
            continue
        for key, value in step_section.items():
            if key not in _STEP_DEFAULT_KEYS:
                continue
            if key == "kwargs" and isinstance(value, Mapping):
                current = merged.setdefault("kwargs", {})
                if isinstance(current, Mapping):
                    merged_kwargs = dict(current)
                else:  # pragma: no cover - defensive branch
                    merged_kwargs = {}
                merged_kwargs.update(dict(value))
                merged["kwargs"] = merged_kwargs
                continue
            if key == "gates" and isinstance(value, Sequence):
                merged[key] = list(value)
                continue
            merged[key] = value
    if "gates" in merged and not isinstance(merged["gates"], list):
        merged["gates"] = []
    return merged


def _index_steps(steps: Sequence[Mapping[str, Any]]) -> tuple[list[str], dict[str, dict[str, Any]]]:
    """Return the step order and index for *steps* ensuring unique identifiers."""

    order: list[str] = []
    index: dict[str, dict[str, Any]] = {}
    for raw in steps:
        if not isinstance(raw, Mapping):
            raise PostprocessingConfigError("step entries must be mappings")
        identifier = raw.get("id")
        if not isinstance(identifier, str) or not identifier:
            raise PostprocessingConfigError("each step must declare a non-empty 'id'")
        if identifier in index:
            raise PostprocessingConfigError(
                f"duplicate step identifier encountered: {identifier}"
            )
        order.append(identifier)
        index[identifier] = dict(raw)
    return order, index


def _merge_mapping(a: Mapping[str, Any], b: Mapping[str, Any]) -> dict[str, Any]:
    result: dict[str, Any] = dict(a)
    for key, value in b.items():
        if key in result and isinstance(result[key], Mapping) and isinstance(value, Mapping):
            result[key] = _merge_mapping(result[key], value)
            continue
        result[key] = value
    return result


def _merge_steps(
    base_steps: Sequence[Mapping[str, Any]],
    override_steps: Sequence[Mapping[str, Any]],
    step_defaults: Mapping[str, Any],
) -> list[dict[str, Any]]:
    base_order, base_index = _index_steps(base_steps)
    override_order, override_index = _index_steps(override_steps)

    merged_steps: list[dict[str, Any]] = []
    for identifier in base_order:
        base = base_index[identifier]
        override = override_index.get(identifier)
        merged = dict(step_defaults)
        merged.update(base)
        if override is not None:
            merged = _merge_mapping(merged, override)
        merged_steps.append(_normalise_step(merged))

    for identifier in override_order:
        if identifier in base_index:
            continue
        override = override_index[identifier]
        merged = dict(step_defaults)
        merged.update(override)
        merged_steps.append(_normalise_step(merged))

    return merged_steps


def _normalise_step(step: Mapping[str, Any]) -> dict[str, Any]:
    """Normalise gate and kwargs entries for *step*."""

    normalised = dict(step)
    gates = normalised.get("gates", [])
    if isinstance(gates, Sequence) and not isinstance(gates, (str, bytes)):
        seen: set[str] = set()
        ordered: list[str] = []
        for gate in gates:
            if isinstance(gate, str) and gate and gate not in seen:
                seen.add(gate)
                ordered.append(gate)
        normalised["gates"] = ordered
    else:
        normalised.pop("gates", None)

    kwargs = normalised.get("kwargs", {})
    if isinstance(kwargs, Mapping):
        normalised["kwargs"] = dict(kwargs)
    else:
        normalised.pop("kwargs", None)
    return normalised


def _ensure_step_callable(step: Mapping[str, Any]) -> str:
    callable_ref = step.get("callable")
    if not isinstance(callable_ref, str) or not callable_ref.strip():
        identifier = step.get("id", "<unknown>")
        raise PostprocessingConfigError(
            f"pipeline step '{identifier}' must define a callable"
        )
    return callable_ref


def _build_steps(
    steps: Sequence[Mapping[str, Any]],
    *,
    flags: Mapping[str, bool],
) -> tuple[PostprocessingStep, ...]:
    built: list[PostprocessingStep] = []
    for raw in steps:
        identifier = raw.get("id")
        if not isinstance(identifier, str):
            raise PostprocessingConfigError("pipeline step is missing an identifier")
        callable_ref = _ensure_step_callable(raw)
        enabled = bool(raw.get("enabled", True))
        continue_on_error = bool(raw.get("continue_on_error", False))
        description = raw.get("description")
        gates_raw = raw.get("gates", [])
        gates: list[str] = []
        if isinstance(gates_raw, Sequence) and not isinstance(gates_raw, (str, bytes)):
            for gate in gates_raw:
                if not isinstance(gate, str) or not gate:
                    continue
                if gate not in flags:
                    raise PostprocessingConfigError(
                        f"pipeline step '{identifier}' references unknown gate '{gate}'"
                    )
                if gate not in gates:
                    gates.append(gate)
        kwargs = raw.get("kwargs", {})
        mapping_kwargs: Mapping[str, Any]
        if isinstance(kwargs, Mapping):
            mapping_kwargs = dict(kwargs)
        else:
            mapping_kwargs = {}
        built.append(
            PostprocessingStep(
                identifier=identifier,
                callable=callable_ref,
                enabled=enabled,
                continue_on_error=continue_on_error,
                description=description if isinstance(description, str) else None,
                gates=tuple(gates),
                kwargs=mapping_kwargs,
            )
        )
    return tuple(built)


def _merge_flags(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, bool]:
    merged: dict[str, bool] = {}
    for source in (base.get("flags"), override.get("flags")):
        if not isinstance(source, Mapping):
            continue
        for key, value in source.items():
            if isinstance(key, str) and key:
                merged[key] = bool(value)
    return merged


def load_postprocessing_pipeline(
    table: str | None,
    *,
    base_dir: Path | str | None = None,
    schema_path: Path | str | None = None,
) -> PostprocessingPipeline:
    """Return the merged pipeline configuration for *table*.

    When *table* is :data:`None`, only the default configuration is returned.
    """

    config_dir = Path(base_dir) if base_dir is not None else POSTPROCESSING_CONFIG_DIR
    schema_file = Path(schema_path) if schema_path is not None else CONFIG_SCHEMA_PATH

    default_path = config_dir / "defaults.yaml"
    if not default_path.exists():
        raise PostprocessingConfigError(
            f"missing default pipeline configuration at {default_path!s}"
        )
    defaults = _load_yaml_mapping(default_path)
    _validate_schema(defaults, schema_file)

    if table is None:
        effective = defaults
    else:
        override_path = config_dir / f"{table}.yaml"
        if not override_path.exists():
            raise PostprocessingConfigError(
                f"missing pipeline configuration for '{table}' at {override_path!s}"
            )
        overrides = _load_yaml_mapping(override_path)
        _validate_schema(overrides, schema_file)
        effective = {
            "version": overrides.get("version", defaults.get("version", 1)),
        }
        effective["flags"] = _merge_flags(defaults, overrides)
        step_defaults = _merge_step_defaults(defaults, overrides)
        base_steps = defaults.get("steps", [])
        override_steps = overrides.get("steps", [])
        if not isinstance(base_steps, Sequence) or isinstance(base_steps, (str, bytes)):
            raise PostprocessingConfigError("default pipeline 'steps' must be a sequence")
        if not isinstance(override_steps, Sequence) or isinstance(
            override_steps, (str, bytes)
        ):
            raise PostprocessingConfigError("override pipeline 'steps' must be a sequence")
        merged_steps = _merge_steps(base_steps, override_steps, step_defaults)
        effective["steps"] = merged_steps
        if step_defaults:
            effective["defaults"] = {"step": step_defaults}
    if table is None:
        step_defaults = _merge_step_defaults(defaults, {})
        flags = _merge_flags(defaults, {})
        steps_raw = defaults.get("steps", [])
        if not isinstance(steps_raw, Sequence) or isinstance(steps_raw, (str, bytes)):
            raise PostprocessingConfigError("default pipeline 'steps' must be a sequence")
        merged_steps = _merge_steps(steps_raw, [], step_defaults)
        effective = {
            "version": defaults.get("version", 1),
            "flags": flags,
            "steps": merged_steps,
        }
        if step_defaults:
            effective["defaults"] = {"step": step_defaults}

    _validate_schema(effective, schema_file)

    flags = effective.get("flags", {})
    if not isinstance(flags, Mapping):
        flags = {}
    steps_raw = effective.get("steps", [])
    if not isinstance(steps_raw, Sequence) or isinstance(steps_raw, (str, bytes)):
        raise PostprocessingConfigError("pipeline 'steps' must be a sequence")
    step_defaults = effective.get("defaults", {}).get("step", {})
    if not isinstance(step_defaults, Mapping):
        step_defaults = {}

    steps = _build_steps(steps_raw, flags=flags)
    version = effective.get("version", 1)
    if not isinstance(version, int):
        raise PostprocessingConfigError("pipeline 'version' must be an integer")

    return PostprocessingPipeline(
        version=version,
        flags=dict(flags),
        steps=steps,
        step_defaults=dict(step_defaults),
    )


def list_postprocessing_tables(
    base_dir: Path | str | None = None,
) -> list[str]:
    """Return available table configuration names in *base_dir*."""

    config_dir = Path(base_dir) if base_dir is not None else POSTPROCESSING_CONFIG_DIR
    if not config_dir.exists():
        return []
    tables: list[str] = []
    for path in sorted(config_dir.glob("*.yaml")):
        if path.name == "defaults.yaml":
            continue
        tables.append(path.stem)
    return tables


def validate_postprocessing_config(
    *,
    tables: Sequence[str] | None = None,
    base_dir: Path | str | None = None,
    schema_path: Path | str | None = None,
) -> list[str]:
    """Validate pipeline configuration files returning processed table names."""

    config_dir = Path(base_dir) if base_dir is not None else POSTPROCESSING_CONFIG_DIR
    schema_file = Path(schema_path) if schema_path is not None else CONFIG_SCHEMA_PATH

    available = list_postprocessing_tables(config_dir)
    target_tables: Sequence[str]
    if tables:
        target_tables = tables
    else:
        target_tables = available

    processed: list[str] = []
    if not target_tables:
        # Still validate defaults so configuration errors are reported early.
        load_postprocessing_pipeline(None, base_dir=config_dir, schema_path=schema_file)
        return processed

    for table in target_tables:
        load_postprocessing_pipeline(table, base_dir=config_dir, schema_path=schema_file)
        processed.append(table)
    return processed


def diff_postprocessing_pipeline(
    table: str,
    *,
    base_dir: Path | str | None = None,
    schema_path: Path | str | None = None,
) -> str:
    """Return a textual diff between defaults and the *table* configuration."""

    config_dir = Path(base_dir) if base_dir is not None else POSTPROCESSING_CONFIG_DIR
    schema_file = Path(schema_path) if schema_path is not None else CONFIG_SCHEMA_PATH

    default_pipeline = load_postprocessing_pipeline(
        None, base_dir=config_dir, schema_path=schema_file
    )
    table_pipeline = load_postprocessing_pipeline(
        table, base_dir=config_dir, schema_path=schema_file
    )

    default_yaml = yaml.safe_dump(
        default_pipeline.as_dict(), sort_keys=False, allow_unicode=True
    ).splitlines(keepends=True)
    table_yaml = yaml.safe_dump(
        table_pipeline.as_dict(), sort_keys=False, allow_unicode=True
    ).splitlines(keepends=True)

    from difflib import unified_diff

    diff = unified_diff(
        default_yaml,
        table_yaml,
        fromfile="defaults",
        tofile=table,
    )
    return "".join(diff)


def _build_parser() -> "argparse.ArgumentParser":  # pragma: no cover - CLI glue
    import argparse

    parser = argparse.ArgumentParser(
        prog="postprocessing-config",
        description="Validate and inspect post-processing pipeline configuration.",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=POSTPROCESSING_CONFIG_DIR,
        help="Directory containing post-processing pipeline YAML files.",
    )
    parser.add_argument(
        "--schema",
        type=Path,
        default=CONFIG_SCHEMA_PATH,
        help="Path to config.schema.json containing the post-processing schema.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_parser = subparsers.add_parser(
        "validate", help="Validate one or more table configurations."
    )
    validate_parser.add_argument(
        "tables",
        nargs="*",
        help="Table names to validate; defaults to all available tables.",
    )

    diff_parser = subparsers.add_parser(
        "diff", help="Show merged configuration diff for a table."
    )
    diff_parser.add_argument("table", help="Table name to diff against defaults.")

    return parser


def main(argv: Sequence[str] | None = None) -> int:  # pragma: no cover - CLI glue
    parser = _build_parser()
    args = parser.parse_args(argv)

    config_dir = args.config_dir
    schema_path = args.schema

    if args.command == "validate":
        tables = args.tables or None
        processed = validate_postprocessing_config(
            tables=tables,
            base_dir=config_dir,
            schema_path=schema_path,
        )
        if processed:
            joined = ", ".join(processed)
            print(f"Validated post-processing configs: {joined}")
        else:
            print("Validated post-processing defaults")
        return 0

    if args.command == "diff":
        diff_text = diff_postprocessing_pipeline(
            args.table, base_dir=config_dir, schema_path=schema_path
        )
        if not diff_text:
            print("No differences between defaults and table configuration")
        else:
            sys.stdout.write(diff_text)
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI glue
    raise SystemExit(main())
