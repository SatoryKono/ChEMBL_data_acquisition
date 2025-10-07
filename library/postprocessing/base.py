"""Shared infrastructure for declarative post-processing pipelines.

This module centralises three concerns required by the refactored helpers:

* ``PipelineContext`` – lightweight container propagated across steps so they
  can access shared state, configuration values and the structured logger.
* ``PostprocessingStep`` – base class enforcing the DataFrame in/out contract
  and exposing the parsed parameters passed via YAML definitions.
* ``PipelineRunner`` – orchestrator that reads declarative pipeline manifests,
  instantiates the registered steps and executes them sequentially while
  emitting deterministic structured logs.

The implementation mirrors the behaviour of the existing CLI orchestrators:
steps are addressed by a table/step name pair, execution is logged with
``postprocess_*`` events, and failures raise rich exceptions so automated
pipelines can surface actionable diagnostics.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Callable, Iterable, Mapping, MutableMapping, Sequence

import pandas as pd
import yaml

from library.common.log import logger as _logger
from .helpers import read_csv_with_fallbacks


DataFrame = pd.DataFrame


class PostprocessingError(RuntimeError):
    """Base exception for post-processing pipeline errors."""


class PipelineConfigurationError(PostprocessingError):
    """Raised when the YAML definition is malformed."""


class StepNotRegisteredError(PostprocessingError):
    """Raised when a pipeline references an unknown step."""


class StepExecutionError(PostprocessingError):
    """Raised when a step fails while transforming the DataFrame."""


def _ensure_dataframe(result: Any, *, table: str, step: str) -> DataFrame:
    if isinstance(result, pd.DataFrame):
        return result
    raise StepExecutionError(
        f"post-processing step '{table}.{step}' did not return a DataFrame"
    )


@dataclass
class PipelineContext:
    """Shared execution context propagated across post-processing steps."""

    table: str
    step: str
    input_path: Path
    output_path: Path | None = None
    params: Mapping[str, Any] = field(default_factory=dict)
    state: MutableMapping[str, Any] | None = None
    logger: Any = _logger

    def __post_init__(self) -> None:
        if self.state is None:
            self.state = {}
        if not isinstance(self.input_path, Path):
            self.input_path = Path(self.input_path)
        if self.output_path is not None and not isinstance(self.output_path, Path):
            self.output_path = Path(self.output_path)

    # ------------------------------------------------------------------
    # Context helpers
    # ------------------------------------------------------------------
    def child(
        self,
        step: str,
        *,
        params: Mapping[str, Any] | None = None,
    ) -> "PipelineContext":
        """Return a copy of the context for ``step`` with merged parameters."""

        merged: dict[str, Any] = dict(self.params)
        if params:
            merged.update(params)
        return PipelineContext(
            table=self.table,
            step=step,
            input_path=self.input_path,
            output_path=self.output_path,
            params=merged,
            state=self.state,
            logger=self.logger,
        )

    # ------------------------------------------------------------------
    # Logging helpers
    # ------------------------------------------------------------------
    def _log(self, level: str, event: str, **kv: Any) -> None:
        payload = {"table": self.table, "step": self.step}
        payload.update(kv)
        self.logger.log(level, event, **payload)

    def info(self, event: str, **kv: Any) -> None:
        self._log("INFO", event, **kv)

    def warn(self, event: str, **kv: Any) -> None:
        self._log("WARN", event, **kv)

    warning = warn

    def error(self, event: str, **kv: Any) -> None:
        self._log("ERROR", event, **kv)


class PostprocessingStep:
    """Base class for DataFrame-to-DataFrame post-processing steps."""

    table_name: str | None = None
    step_name: str | None = None
    description: str | None = None

    def __init__(self, *, params: Mapping[str, Any] | None = None) -> None:
        self._params = MappingProxyType(dict(params or {}))

    # ------------------------------------------------------------------
    # Metadata helpers
    # ------------------------------------------------------------------
    @property
    def name(self) -> str:
        return self.step_name or self.__class__.__name__

    @property
    def params(self) -> Mapping[str, Any]:
        return self._params

    # ------------------------------------------------------------------
    # Execution entry point
    # ------------------------------------------------------------------
    def __call__(self, frame: DataFrame, context: PipelineContext) -> DataFrame:
        result = self.transform(frame, context)
        return _ensure_dataframe(result, table=context.table, step=self.name)

    def transform(self, frame: DataFrame, context: PipelineContext) -> DataFrame:
        """Return a transformed DataFrame.

        Sub-classes must override this method. The base implementation raises
        :class:`NotImplementedError` to highlight missing overrides early during
        development.
        """

        raise NotImplementedError("PostprocessingStep.transform must be overridden")


_STEP_REGISTRY: dict[str, dict[str, type[PostprocessingStep]]] = {}


def register_step(
    *, table: str, name: str | None = None, description: str | None = None
) -> Callable[[type[PostprocessingStep]], type[PostprocessingStep]]:
    """Register a :class:`PostprocessingStep` for ``table``.

    The decorator enforces uniqueness of the ``table``/``name`` pair and augments
    the class with the resolved metadata so callers can introspect the registry.
    """

    def decorator(cls: type[PostprocessingStep]) -> type[PostprocessingStep]:
        if not issubclass(cls, PostprocessingStep):
            raise TypeError("registered object must subclass PostprocessingStep")

        step_name = name or getattr(cls, "step_name", None) or cls.__name__
        table_registry = _STEP_REGISTRY.setdefault(table, {})
        if step_name in table_registry:
            raise ValueError(f"post-processing step already registered: {table}.{step_name}")

        cls.table_name = table
        cls.step_name = step_name
        cls.description = description
        table_registry[step_name] = cls
        return cls

    return decorator


def get_registered_step(table: str, name: str) -> type[PostprocessingStep]:
    """Return the registered class for ``table``/``name``."""

    table_registry = _STEP_REGISTRY.get(table)
    if table_registry is None or name not in table_registry:
        raise StepNotRegisteredError(f"post-processing step not found: {table}.{name}")
    return table_registry[name]


def iter_registered_steps(table: str | None = None) -> Iterable[PostprocessingStep]:
    """Yield instantiated steps for ``table`` or all tables if ``None``."""

    tables: Iterable[str]
    if table is None:
        tables = _STEP_REGISTRY.keys()
    else:
        tables = (table,)
    for table_name in tables:
        for step_cls in _STEP_REGISTRY.get(table_name, {}).values():
            yield step_cls()


@dataclass(frozen=True)
class _StepDefinition:
    name: str
    params: Mapping[str, Any]


@dataclass(frozen=True)
class _PipelineDefinition:
    table: str
    input_path: Path
    steps: tuple[_StepDefinition, ...]
    output_path: Path | None = None
    params: Mapping[str, Any] = field(default_factory=dict)


class PipelineRunner:
    """Execute registered post-processing steps declared in YAML manifests."""

    def __init__(
        self,
        pipelines: Sequence[_PipelineDefinition],
        *,
        registry: Mapping[str, Mapping[str, type[PostprocessingStep]]] | None = None,
    ) -> None:
        self._pipelines = tuple(pipelines)
        self._registry = registry or _STEP_REGISTRY

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------
    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineRunner":
        manifest_path = Path(path)
        data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        if not data:
            raise PipelineConfigurationError("post-processing manifest is empty")
        pipelines_raw = data.get("pipelines")
        if not pipelines_raw:
            raise PipelineConfigurationError("manifest missing 'pipelines' section")

        base_dir = manifest_path.parent
        if isinstance(pipelines_raw, Mapping):
            entries: list[Mapping[str, Any]] = []
            for table_name, cfg in pipelines_raw.items():
                if not isinstance(cfg, Mapping):
                    raise PipelineConfigurationError(
                        "pipeline configuration entries must be mappings"
                    )
                merged: dict[str, Any] = {"table": table_name}
                merged.update(cfg)
                entries.append(merged)
        elif isinstance(pipelines_raw, Sequence) and not isinstance(
            pipelines_raw, (str, bytes)
        ):
            entries = list(pipelines_raw)  # type: ignore[arg-type]
        else:
            raise PipelineConfigurationError(
                "'pipelines' must be a sequence or mapping of pipeline definitions"
            )

        definitions = tuple(
            cls._coerce_pipeline(entry, base_dir=base_dir) for entry in entries
        )
        return cls(definitions)

    @staticmethod
    def _coerce_pipeline(
        entry: Mapping[str, Any],
        *,
        base_dir: Path,
    ) -> _PipelineDefinition:
        if not isinstance(entry, Mapping):
            raise PipelineConfigurationError("each pipeline entry must be a mapping")

        table = entry.get("table") or entry.get("name")
        if not table:
            raise PipelineConfigurationError("pipeline entry missing 'table'")

        input_value = entry.get("input")
        if not input_value:
            raise PipelineConfigurationError(f"pipeline '{table}' missing 'input'")

        input_path = Path(input_value)
        if not input_path.is_absolute():
            input_path = (base_dir / input_path).resolve()

        output_value = entry.get("output")
        output_path: Path | None
        if output_value:
            output_path = Path(output_value)
            if not output_path.is_absolute():
                output_path = (base_dir / output_path).resolve()
        else:
            output_path = None

        params = entry.get("params") or {}
        if not isinstance(params, Mapping):
            raise PipelineConfigurationError(
                f"pipeline '{table}' field 'params' must be a mapping"
            )

        steps_raw = entry.get("steps")
        if not steps_raw:
            raise PipelineConfigurationError(f"pipeline '{table}' missing 'steps'")

        steps: list[_StepDefinition] = []
        for item in steps_raw:
            if isinstance(item, str):
                step_name = item
                step_params: Mapping[str, Any] = {}
            elif isinstance(item, Mapping):
                step_name = item.get("name") or item.get("step")
                if not step_name:
                    raise PipelineConfigurationError(
                        f"pipeline '{table}' step missing 'name'"
                    )
                declared_params = item.get("params") or item.get("options") or {}
                if not isinstance(declared_params, Mapping):
                    raise PipelineConfigurationError(
                        f"pipeline '{table}' step '{step_name}' params must be a mapping"
                    )
                step_params = declared_params
            else:
                raise PipelineConfigurationError(
                    f"pipeline '{table}' contains invalid step definition"
                )
            steps.append(_StepDefinition(name=step_name, params=step_params))

        return _PipelineDefinition(
            table=str(table),
            input_path=input_path,
            output_path=output_path,
            params=dict(params),
            steps=tuple(steps),
        )

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def run(self) -> dict[str, DataFrame]:
        results: dict[str, DataFrame] = {}
        for definition in self._pipelines:
            pipeline_logger = _logger.bind(table=definition.table)
            pipeline_logger.info(
                "postprocess_pipeline_start", input=str(definition.input_path)
            )

            if definition.table not in self._registry:
                raise StepNotRegisteredError(
                    f"no steps registered for table '{definition.table}'"
                )

            frame = read_csv_with_fallbacks(definition.input_path)

            base_context = PipelineContext(
                table=definition.table,
                step="pipeline",
                input_path=definition.input_path,
                output_path=definition.output_path,
                params=definition.params,
                state={},
                logger=pipeline_logger,
            )

            for step_definition in definition.steps:
                step_cls = self._resolve_step(definition.table, step_definition.name)
                step_context = base_context.child(
                    step_cls.step_name or step_definition.name,
                    params=step_definition.params,
                )
                step = step_cls(params=step_context.params)

                step_context.info(
                    "postprocess_step_start", rows=len(frame), input=str(definition.input_path)
                )
                try:
                    frame = step(frame, step_context)
                except Exception as exc:  # pragma: no cover - defensive wrapper
                    step_context.error(
                        "postprocess_step_failed",
                        error=str(exc),
                        exception=exc.__class__.__name__,
                    )
                    raise StepExecutionError(
                        f"post-processing step '{definition.table}.{step.name}' failed"
                    ) from exc
                else:
                    step_context.info(
                        "postprocess_step_done",
                        rows=len(frame),
                        output=str(definition.output_path) if definition.output_path else None,
                    )

            pipeline_logger.info(
                "postprocess_pipeline_done",
                rows=len(frame),
                output=str(definition.output_path) if definition.output_path else None,
            )
            results[definition.table] = frame
        return results

    def _resolve_step(
        self, table: str, name: str
    ) -> type[PostprocessingStep]:  # pragma: no cover - thin wrapper
        table_registry = self._registry.get(table)
        if table_registry is None or name not in table_registry:
            raise StepNotRegisteredError(f"post-processing step not found: {table}.{name}")
        return table_registry[name]


__all__ = [
    "DataFrame",
    "PipelineContext",
    "PipelineRunner",
    "PipelineConfigurationError",
    "PostprocessingError",
    "PostprocessingStep",
    "StepExecutionError",
    "StepNotRegisteredError",
    "get_registered_step",
    "iter_registered_steps",
    "register_step",
]

