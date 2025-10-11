from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any, TypedDict, cast

import yaml


class PipelineStepFlags(TypedDict, total=False):
    """Optional capability flags supported by a pipeline step."""

    dry_run: bool


class PipelineStepDefinition(TypedDict, total=False):
    """Declarative description of a pipeline stage."""

    name: str
    callable: str
    input: str
    output: str
    subcommand: str | None
    output_flag: str
    extra_args: list[str]
    flags: PipelineStepFlags
    depends_on: list[str]
    produces: list[str]
    consumes: list[str]


@dataclass(frozen=True)
class PipelineStep:
    """Runtime representation of a pipeline stage."""

    name: str
    main: Callable[[Sequence[str] | None], int]
    input_filename: str
    output_stem: str
    subcommand: str | None = None
    output_flag: str = "--final-out"
    extra_args: tuple[str, ...] = ()
    supports_dry_run: bool = False
    depends_on: tuple[str, ...] = ()
    produces: tuple[str, ...] = ()
    consumes: tuple[str, ...] = ()

    def build_arguments(self, cfg: Any, output_path: Path | None = None) -> list[str]:
        """Return CLI arguments forwarded to the wrapped ``main`` function."""

        input_csv = cfg.input_path(self.name)
        output_csv = (
            output_path if output_path is not None else cfg.output_path(self.name)
        )
        args = ["--config", str(cfg.config_path), "--input", str(input_csv)]
        args.extend([self.output_flag, str(output_csv)])
        args.extend(["--log-level", cfg.log_level])
        limit = getattr(cfg, "limit", None)
        if limit is not None:
            args.extend(["--limit", str(limit)])
        if getattr(cfg, "force", False):
            args.append("--force")
        if getattr(cfg, "skip_existing", False):
            args.append("--skip-existing")
        if getattr(cfg, "dry_run", False) and self.supports_dry_run:
            args.append("--dry-run")
        if self.extra_args:
            args = [*self.extra_args, *args]
        if self.subcommand is not None:
            return [self.subcommand, *args]
        return args

    def expected_output(self, cfg: Any) -> Path:
        """Return the path where the pipeline will create its CSV artefact."""

        output = cfg.output_path(self.name)
        return output if isinstance(output, Path) else Path(output)

    def required_input(self, cfg: Any) -> Path:
        """Return the CSV that the pipeline expects as input."""

        input_path = cfg.input_path(self.name)
        return input_path if isinstance(input_path, Path) else Path(input_path)


_DEFAULT_DEFINITIONS: tuple[PipelineStepDefinition, ...] = (
    {
        "name": "document",
        "callable": "library.cli.commands.get_document_data:main",
        "input": "document.csv",
        "output": "documents",
        "extra_args": ["--mode", "all"],
        "produces": ["documents"],
    },
    {
        "name": "target",
        "callable": "library.cli.commands.get_target_data:main",
        "input": "target.csv",
        "output": "targets",
        "subcommand": "all",
        "output_flag": "--final-out",
        "produces": ["targets"],
    },
    {
        "name": "assay",
        "callable": "library.cli.commands.get_assay_data:main",
        "input": "assay.csv",
        "output": "assays",
        "produces": ["assays"],
    },
    {
        "name": "testitem",
        "callable": "library.cli.commands.get_testitem_data:main",
        "input": "testitem.csv",
        "output": "testitem",
        "produces": ["testitem"],
    },
    {
        "name": "activity",
        "callable": "library.cli.entrypoints.activity:main",
        "input": "activity.csv",
        "output": "activity",
        "flags": {"dry_run": True},
        "produces": ["activity", "activities"],
        "consumes": ["documents", "targets", "assays", "testitem"],
    },
)


def load_pipeline_registry(
    source: (
        Path | str | Iterable[PipelineStepDefinition] | Mapping[str, object] | None
    ) = None,
) -> tuple[PipelineStep, ...]:
    """Return pipeline steps loaded from ``source`` or the default registry."""

    definitions: Iterable[PipelineStepDefinition]
    if source is None:
        definitions = _DEFAULT_DEFINITIONS
    elif isinstance(source, str | Path):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"registry file not found: {path}")
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        definitions = _coerce_definitions(data)
    elif isinstance(source, Mapping):
        definitions = _coerce_definitions(source)
    else:
        definitions = source
    return tuple(_build_step(entry) for entry in definitions)


def _coerce_definitions(data: object) -> Iterable[PipelineStepDefinition]:
    if data is None:
        return ()
    if isinstance(data, Mapping):
        pipelines = data.get("pipelines")
        if pipelines is None:
            raise ValueError("registry mapping must contain a 'pipelines' key")
        data = pipelines
    if not isinstance(data, Iterable) or isinstance(data, str | bytes):
        raise TypeError("registry definitions must be an iterable of mappings")
    definitions: list[PipelineStepDefinition] = []
    for entry in data:
        if not isinstance(entry, Mapping):
            raise TypeError("each pipeline definition must be a mapping")
        definitions.append(cast(PipelineStepDefinition, dict(entry)))
    return definitions


def _coerce_string_sequence(value: object, *, field: str) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    if isinstance(value, Iterable) and not isinstance(value, bytes | bytearray):
        result: list[str] = []
        for item in value:
            if not isinstance(item, str):
                raise TypeError(
                    f"pipeline definition field '{field}' must contain strings only"
                )
            result.append(item)
        return tuple(result)
    raise TypeError(
        f"pipeline definition field '{field}' must be a string or an iterable of strings"
    )


def _build_step(entry: PipelineStepDefinition) -> PipelineStep:
    name = entry.get("name")
    if not name:
        raise ValueError("pipeline definition missing 'name'")
    dotted = entry.get("callable")
    if not dotted:
        raise ValueError(f"pipeline '{name}' missing 'callable'")
    main = _resolve_callable(dotted)
    input_filename = entry.get("input") or f"{name}.csv"
    output_stem = entry.get("output") or name
    subcommand = entry.get("subcommand")
    output_flag = entry.get("output_flag") or "--final-out"
    extra_args = tuple(entry.get("extra_args", ()))
    flags = entry.get("flags", {})
    supports_dry_run = bool(flags.get("dry_run", False))
    depends_on = _coerce_string_sequence(entry.get("depends_on"), field="depends_on")
    produces_raw = entry.get("produces", (output_stem,))
    produces = _coerce_string_sequence(produces_raw, field="produces")
    consumes = _coerce_string_sequence(entry.get("consumes"), field="consumes")
    return PipelineStep(
        name=name,
        main=main,
        input_filename=input_filename,
        output_stem=output_stem,
        subcommand=subcommand,
        output_flag=output_flag,
        extra_args=extra_args,
        supports_dry_run=supports_dry_run,
        depends_on=depends_on,
        produces=produces,
        consumes=consumes,
    )


def _resolve_callable(dotted: str) -> Callable[[Sequence[str] | None], int]:
    module_path, _, attribute = dotted.partition(":")
    if not module_path or not attribute:
        raise ValueError(f"invalid callable reference: {dotted}")
    module = import_module(module_path)
    target = module
    for part in attribute.split("."):
        if not hasattr(target, part):
            raise AttributeError(f"callable '{dotted}' not found")
        target = getattr(target, part)
    if not callable(target):
        raise TypeError(f"resolved object for '{dotted}' is not callable")
    return target


__all__ = [
    "PipelineStep",
    "PipelineStepDefinition",
    "PipelineStepFlags",
    "load_pipeline_registry",
]
