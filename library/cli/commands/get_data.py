"""Orchestrate all ChEMBL data acquisition pipelines via a unified CLI interface.

This module provides a single entry point that wires the existing acquisition
scripts together. ``get_data.py`` accepts common options for locating input and
output artefacts, resolves those paths once, and then invokes the following
modules in order while forwarding the computed arguments:

1. :mod:`scripts.get_document_data`
2. :mod:`scripts.get_target_data`
3. :mod:`scripts.get_assay_data`
4. :mod:`scripts.get_testitem_data`
5. :mod:`scripts.get_activity_data`

Each pipeline receives consistent ``--config``, ``--input`` and ``--final-out``
values alongside shared logging options. The command line interface exposes
parameters for the base data directory, distinct input/output sub-directories,
a date prefix for generated filenames and optional execution flags. The helper
functions defined below encapsulate the argument preparation for every step so
that the pipelines can be executed programmatically from other callers as well.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import time
from collections import deque
from collections.abc import Callable, Iterable, ItemsView, KeysView, Mapping, Sequence, ValuesView
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from fnmatch import fnmatch
from heapq import heappop, heappush
from pathlib import Path
from typing import IO, Any, Optional, TypeVar, Generic, Iterator

from pydantic import BaseModel, ConfigDict, ValidationError

from library.cli.logging import setup_cli_logging
from library.cli_utils import resolve_invocation
from library.clients import ChemblClient
from library.common.logging_setup import Logger, LoggerConfig, configure_logger
from library.config import (
    DEFAULT_CONFIG_PATH,
    Config,
    ConfigError,
    ConfigLoaderError,
    ensure_dirs,
    load_config,
    print_config,
)
from library.orchestration import ETLContext
from library.orchestration.workflow import (
    PreparedPipelineStep,
    StepExecutionResult,
    execute_workflow,
    temporary_output_path,
)
from library.pipelines.activity import (
    ActivityPipelineOptions,
)
from library.pipelines.activity import (
    run_pipeline as run_activity_pipeline,
)
from library.pipelines.assay import (
    AssayPipelineOptions,
)
from library.pipelines.assay import (
    run_pipeline as run_assay_pipeline,
)
from library.pipelines.common import PipelineRunResult
from library.pipelines.document import (
    DocumentPipelineOptions,
)
from library.pipelines.document import (
    run_pipeline as run_document_pipeline,
)
from library.pipelines.registry import PipelineStep, load_pipeline_registry
from library.pipelines.target import (
    TargetPipelineOptions,
)
from library.pipelines.target import (
    run_pipeline as run_target_pipeline,
)
from library.pipelines.testitem import (
    TestitemPipelineOptions,
)
from library.pipelines.testitem import (
    run_pipeline as run_testitem_pipeline,
)
from library.postprocess.common import (
    SUPPORTED_TABLES as POSTPROCESS_SUPPORTED_TABLES,
)
from library.postprocess.common import (
    PostprocessingPipelineConfig,
    PostprocessResult,
    run_postprocessing_pipeline,
)
from library.postprocess.common import (
    get_csv_runtime_config as get_postprocess_csv_config,
)
from library.postprocess.common import (
    get_pipeline_config as load_postprocess_pipeline_config,
)
from library.postprocessing.activities import (
    ACTIVITY_SCHEMA,
    validate_activities,
)
from library.postprocessing.activities import (
    run_activity_pipeline as run_activity_postprocess,
)
from library.postprocessing.assays import (
    ASSAY_SCHEMA,
    validate_assays,
)
from library.postprocessing.assays import (
    run_assay_pipeline as run_assay_postprocess,
)
from library.postprocessing.documents import (
    DOCUMENT_SCHEMA,
    validate_documents,
)
from library.postprocessing.documents import (
    run_document_pipeline as run_document_postprocess,
)
from library.postprocessing.targets import (
    TARGET_SCHEMA,
    validate_targets,
)
from library.postprocessing.targets import (
    run_target_pipeline as run_target_postprocess,
)
from library.postprocessing.testitem import (
    TESTITEM_SCHEMA,
    validate_testitems,
)
from library.postprocessing.testitem import (
    run_testitem_pipeline as run_testitem_postprocess,
)
from library.reporting.run_manifest import load_output_report, merge_run_output

_LOGGER: Logger = configure_logger(LoggerConfig())


StepValueT = TypeVar("StepValueT")


class _PipelineMapping(BaseModel, Mapping[str, StepValueT], Generic[StepValueT]):
    """Immutable mapping of pipeline step names to typed values."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    entries: dict[str, StepValueT]

    def __getitem__(self, key: str) -> StepValueT:
        return self.entries[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def __contains__(self, key: object) -> bool:
        return key in self.entries

    def get(self, key: str, default: StepValueT | None = None) -> StepValueT | None:
        return self.entries.get(key, default)

    def items(self) -> ItemsView[str, StepValueT]:
        return self.entries.items()

    def keys(self) -> KeysView[str]:
        return self.entries.keys()

    def values(self) -> ValuesView[StepValueT]:
        return self.entries.values()

    def to_dict(self) -> dict[str, StepValueT]:
        """Return a shallow copy of the underlying mapping."""

        return dict(self.entries)

    @classmethod
    def from_mapping(
        cls, mapping: Mapping[str, StepValueT] | "_PipelineMapping[StepValueT]"
    ) -> "_PipelineMapping[StepValueT]":
        if isinstance(mapping, cls):
            return mapping
        return cls(entries=dict(mapping))


class PipelineInputFiles(_PipelineMapping[str]):
    """Mapping of pipeline step names to their expected input filenames."""


class PipelineOutputStems(_PipelineMapping[str]):
    """Mapping of pipeline step names to the stem used for output artefacts."""


class PipelineSubcommands(_PipelineMapping[str | None]):
    """Mapping of pipeline step names to optional CLI sub-commands."""

    def get(self, key: str, default: str | None = None) -> str | None:  # type: ignore[override]
        return super().get(key, default)


DEFAULT_PIPELINE_STEPS: tuple[PipelineStep, ...] = load_pipeline_registry()
DEFAULT_INPUT_FILES: PipelineInputFiles = PipelineInputFiles.from_mapping(
    {step.name: step.input_filename for step in DEFAULT_PIPELINE_STEPS}
)
DEFAULT_OUTPUT_STEMS: PipelineOutputStems = PipelineOutputStems.from_mapping(
    {step.name: step.output_stem for step in DEFAULT_PIPELINE_STEPS}
)
DEFAULT_SUBCOMMANDS: PipelineSubcommands = PipelineSubcommands.from_mapping(
    {step.name: step.subcommand for step in DEFAULT_PIPELINE_STEPS}
)

# Backwards compatibility for existing callers/tests that patch the legacy names.
_PIPELINE_STEPS = DEFAULT_PIPELINE_STEPS
_DEFAULT_INPUT_FILES = DEFAULT_INPUT_FILES
_DEFAULT_OUTPUT_STEMS = DEFAULT_OUTPUT_STEMS
_DEFAULT_SUBCOMMANDS = DEFAULT_SUBCOMMANDS


@dataclass(frozen=True)
class _PostprocessHandlers:
    runner: Callable[..., tuple[Any, Any]]
    validator: Callable[..., Any]
    schema: Any


_POSTPROCESS_HANDLERS: dict[str, _PostprocessHandlers] = {
    "activities": _PostprocessHandlers(
        runner=run_activity_postprocess,
        validator=validate_activities,
        schema=ACTIVITY_SCHEMA,
    ),
    "assays": _PostprocessHandlers(
        runner=run_assay_postprocess,
        validator=validate_assays,
        schema=ASSAY_SCHEMA,
    ),
    "documents": _PostprocessHandlers(
        runner=run_document_postprocess,
        validator=validate_documents,
        schema=DOCUMENT_SCHEMA,
    ),
    "targets": _PostprocessHandlers(
        runner=run_target_postprocess,
        validator=validate_targets,
        schema=TARGET_SCHEMA,
    ),
    "testitems": _PostprocessHandlers(
        runner=run_testitem_postprocess,
        validator=validate_testitems,
        schema=TESTITEM_SCHEMA,
    ),
}


_UNLINK_MAX_ATTEMPTS = 5
_UNLINK_RETRY_SLEEP_SECONDS = 0.1
_WINDOWS_SHARING_VIOLATION = 32

_DEFAULT_DATE_ENV = "CHEMBL_DA_DEFAULT_DATE"
_DEFAULT_DATE_ENV_ALIAS = "CHEMBL_DA_DEFAULT_DATE_PREFIX"
_DEFAULT_DATE_PREFIX = "19700101"
_RUN_ID_ENV = "CHEMBL_DA_RUN_ID"


OptionsT = TypeVar("OptionsT")


@dataclass(frozen=True)
class PipelineApi(Generic[OptionsT]):
    """Describe how to build options and execute a pipeline programmatically."""

    build_options: Callable[[PipelineRunConfig, Path, Path], OptionsT]
    runner: Callable[[Config, OptionsT], PipelineRunResult]


def _build_document_options(
    cfg: PipelineRunConfig, input_path: Path, output_path: Path
) -> DocumentPipelineOptions:
    mode = cfg.subcommand_for("document") or "all"
    return DocumentPipelineOptions(
        input_csv=input_path,
        output_csv=output_path,
        mode=mode,
        limit=cfg.limit,
        force=cfg.force,
        skip_existing=cfg.skip_existing,
        rerun_postprocess=cfg.rerun_postprocess,
        date_prefix=cfg.date_prefix,
        output_stem=cfg.output_stems.get("document"),
    )


def _build_target_options(
    cfg: PipelineRunConfig, input_path: Path, output_path: Path
) -> TargetPipelineOptions:
    command = cfg.subcommand_for("target") or "all"
    return TargetPipelineOptions(
        input_csv=input_path,
        output_csv=output_path,
        command=command,
        limit=cfg.limit,
        force=cfg.force,
        skip_existing=cfg.skip_existing,
    )


def _build_assay_options(
    cfg: PipelineRunConfig, input_path: Path, output_path: Path
) -> AssayPipelineOptions:
    return AssayPipelineOptions(
        input_csv=input_path,
        output_csv=output_path,
        limit=cfg.limit,
        force=cfg.force,
        skip_existing=cfg.skip_existing,
    )


def _build_testitem_options(
    cfg: PipelineRunConfig, input_path: Path, output_path: Path
) -> TestitemPipelineOptions:
    return TestitemPipelineOptions(
        input_csv=input_path,
        output_csv=output_path,
        limit=cfg.limit,
        offset=0,
    )


def _build_activity_options(
    cfg: PipelineRunConfig, input_path: Path, output_path: Path
) -> ActivityPipelineOptions:
    return ActivityPipelineOptions(
        input_csv=input_path,
        output_csv=output_path,
        limit=cfg.limit,
        force=cfg.force,
        skip_existing=cfg.skip_existing,
        dry_run=cfg.dry_run,
    )


_PIPELINE_APIS: Mapping[str, PipelineApi[Any]] = {
    "document": PipelineApi[DocumentPipelineOptions](
        _build_document_options, run_document_pipeline
    ),
    "target": PipelineApi[TargetPipelineOptions](
        _build_target_options, run_target_pipeline
    ),
    "assay": PipelineApi[AssayPipelineOptions](
        _build_assay_options, run_assay_pipeline
    ),
    "testitem": PipelineApi[TestitemPipelineOptions](
        _build_testitem_options, run_testitem_pipeline
    ),
    "activity": PipelineApi[ActivityPipelineOptions](
        _build_activity_options, run_activity_pipeline
    ),
}


@dataclass(frozen=True)
class PipelineRunConfig:
    """Resolved configuration shared across pipeline steps."""

    base_path: Path
    input_dir: Path
    output_dir: Path
    config_path: Path
    date_prefix: str
    log_level: str
    limit: int | None
    force: bool
    skip_existing: bool
    dry_run: bool
    input_files: PipelineInputFiles
    output_stems: PipelineOutputStems
    subcommands: PipelineSubcommands
    rerun_postprocess: bool = False

    def input_path(self, name: str) -> Path:
        """Return the fully resolved path for ``name`` in the input directory."""

        filename = self.input_files[name]
        return self.input_dir / filename

    def output_path(self, name: str) -> Path:
        """Return the fully resolved path for ``name`` in the output directory."""

        stem = self.output_stems[name]
        stem_path = Path(stem)

        # Allow overrides to provide explicit filenames (e.g. ``output.targets.csv``)
        # or nested locations. When the stem includes a suffix we treat it as a
        # concrete path instead of appending the canonical prefix/suffix.
        if stem_path.suffix:
            if stem_path.is_absolute():
                return stem_path
            return self.output_dir / stem_path

        filename = f"output.{stem}_{self.date_prefix}.csv"
        return self.output_dir / filename

    def subcommand_for(self, name: str) -> str | None:
        """Return the configured subcommand for ``name`` if available."""

        return self.subcommands.get(name)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "input_files", PipelineInputFiles.from_mapping(self.input_files)
        )
        object.__setattr__(
            self, "output_stems", PipelineOutputStems.from_mapping(self.output_stems)
        )
        object.__setattr__(
            self, "subcommands", PipelineSubcommands.from_mapping(self.subcommands)
        )


def _resolve_path(base: Path, candidate: Path) -> Path:
    """Return ``candidate`` if absolute, otherwise relative to ``base``."""

    expanded = candidate.expanduser()
    if expanded.is_absolute():
        return expanded.resolve()
    return (base / expanded).resolve()


def _resolve_consumed_artifact_path(cfg: PipelineRunConfig, artefact: str) -> Path:
    """Return the filesystem path associated with a consumed artefact name."""

    candidate = Path(artefact)
    if candidate.is_absolute():
        return candidate
    if candidate.suffix:
        return cfg.output_dir / candidate
    return cfg.output_dir / f"output.{candidate.name}_{cfg.date_prefix}.csv"


def _parse_overrides(
    values: Sequence[str] | None,
    *,
    allow_empty_value: bool = False,
) -> dict[str, str]:
    """Parse ``STEP=value`` pairs from the CLI into a dictionary."""

    overrides: dict[str, str] = {}
    if not values:
        return overrides
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"invalid override format: {raw!r}")
        name, value = raw.split("=", 1)
        key = name.strip()
        if not key:
            raise ValueError(f"override missing step name: {raw!r}")
        if not value and not allow_empty_value:
            raise ValueError(f"override missing value for step {key!r}")
        overrides[key] = value
    return overrides


def _resolve_pipeline_steps(
    args: argparse.Namespace | None = None,
) -> tuple[PipelineStep, ...]:
    """Load pipeline definitions applying CLI overrides when provided."""

    registry_source = getattr(args, "pipeline_registry", None)
    steps = (
        load_pipeline_registry(registry_source)
        if registry_source is not None
        else DEFAULT_PIPELINE_STEPS
    )
    steps = tuple(steps)

    input_overrides = _parse_overrides(getattr(args, "override_input", None))
    output_overrides = _parse_overrides(getattr(args, "override_output_stem", None))
    subcommand_overrides = _parse_overrides(
        getattr(args, "override_subcommand", None), allow_empty_value=True
    )

    if not input_overrides and not output_overrides and not subcommand_overrides:
        return steps

    known_steps = {step.name for step in steps}
    _validate_override_keys(input_overrides, known_steps, "input")
    _validate_override_keys(output_overrides, known_steps, "output")
    _validate_override_keys(subcommand_overrides, known_steps, "subcommand")

    mutated: list[PipelineStep] = []
    for step in steps:
        updated = step
        if step.name in input_overrides:
            updated = replace(updated, input_filename=input_overrides[step.name])
        if step.name in output_overrides:
            updated = replace(updated, output_stem=output_overrides[step.name])
        if step.name in subcommand_overrides:
            raw_value = subcommand_overrides[step.name]
            new_subcommand = raw_value if raw_value else None
            updated = replace(updated, subcommand=new_subcommand)
        mutated.append(updated)
    return tuple(mutated)


def _validate_override_keys(
    overrides: Mapping[str, str],
    known: set[str],
    kind: str,
) -> None:
    """Ensure override keys reference known steps."""

    unknown = sorted(set(overrides) - known)
    if unknown:
        joined = ", ".join(unknown)
        raise ValueError(f"unknown {kind} override for step(s): {joined}")


def _normalise_date_prefix(value: str) -> str:
    """Validate and normalise a date prefix value."""

    candidate = value.strip()
    if not candidate:
        raise ValueError("date prefix must not be empty when resolved automatically")
    if not candidate.isdigit() or len(candidate) != 8:
        raise ValueError(
            "date prefix must be an eight digit YYYYMMDD string, got" f" {value!r}"
        )
    return candidate


def _resolve_default_date_prefix(
    args: argparse.Namespace,
    *,
    base_path: Path,
) -> str:
    """Return a deterministic default date prefix for CLI invocations."""

    for env_name in (_DEFAULT_DATE_ENV, _DEFAULT_DATE_ENV_ALIAS):
        env_value = os.environ.get(env_name)
        if env_value is not None:
            return _normalise_date_prefix(env_value)

    config_candidate = getattr(args, "config", None) or DEFAULT_CONFIG_PATH
    config_path = Path(config_candidate).expanduser()

    if config_path.exists():
        try:
            base_for_config = base_path if base_path.exists() else None
            config_obj = load_config(config_path, base_path=base_for_config)
        except (ConfigError, ConfigLoaderError, ValidationError, OSError):
            pass
        else:
            local_cfg = getattr(config_obj, "local", None)
            io_cfg = getattr(local_cfg, "io", None)
            default_prefix = getattr(io_cfg, "default_date_prefix", None)
            if isinstance(default_prefix, str) and default_prefix.strip():
                return _normalise_date_prefix(default_prefix)

    return _DEFAULT_DATE_PREFIX


def _ensure_date_prefix(args: argparse.Namespace, *, base_path: Path) -> str:
    """Ensure ``args.date_prefix`` is populated with a deterministic value."""

    current = getattr(args, "date_prefix", None)
    if isinstance(current, str):
        stripped = current.strip()
        if not stripped:
            raise ValueError("--date must not be an empty string")
        args.date_prefix = stripped
        return stripped

    default_prefix = _resolve_default_date_prefix(args, base_path=base_path)
    args.date_prefix = default_prefix
    return default_prefix


def _normalise_descriptor_path(path: Path) -> str:
    try:
        return str(path.expanduser().resolve())
    except OSError:
        return str(path.expanduser())


def _canonical_run_descriptor(args: argparse.Namespace, *, base_path: Path) -> str:
    """Return a canonical description of the invocation for run hashing."""

    invocation = getattr(args, "invocation", None)
    if isinstance(invocation, Sequence) and invocation:
        parts = [str(part) for part in invocation]
    else:
        parts = [str(part) for part in resolve_invocation("get_data", None)]

    resolved_base = base_path.expanduser().resolve()
    parts.append(f"base_path={_normalise_descriptor_path(resolved_base)}")

    input_dir = getattr(args, "input_dir", Path("input"))
    output_dir = getattr(args, "output_dir", Path("output"))
    resolved_input = _resolve_path(resolved_base, Path(input_dir))
    resolved_output = _resolve_path(resolved_base, Path(output_dir))
    parts.append(f"input_dir={_normalise_descriptor_path(resolved_input)}")
    parts.append(f"output_dir={_normalise_descriptor_path(resolved_output)}")

    config_candidate = getattr(args, "config", None) or DEFAULT_CONFIG_PATH
    config_path = Path(config_candidate)
    parts.append(f"config={_normalise_descriptor_path(config_path)}")

    registry_value = getattr(args, "pipeline_registry", None)
    if registry_value not in (None, argparse.SUPPRESS):
        registry_path = Path(registry_value)
        parts.append(f"pipeline_registry={_normalise_descriptor_path(registry_path)}")

    date_prefix = getattr(args, "date_prefix", None)
    if date_prefix is not None:
        parts.append(f"date_prefix={date_prefix}")

    desired_level = getattr(args, "log_level", "INFO")
    parts.append(f"log_level={str(desired_level).upper()}")

    limit_value = getattr(args, "limit", None)
    parts.append(f"limit={limit_value}")

    bool_fields = (
        "force",
        "skip_existing",
        "dry_run",
        "verbose",
        "print_config",
    )
    for field in bool_fields:
        parts.append(f"{field}={bool(getattr(args, field, False))}")

    input_overrides = _parse_overrides(getattr(args, "override_input", None))
    if input_overrides:
        resolved_items: list[str] = []
        for key, value in sorted(input_overrides.items()):
            candidate = Path(value)
            if candidate.expanduser().is_absolute():
                resolved_value = _normalise_descriptor_path(candidate)
            else:
                resolved_value = _normalise_descriptor_path(resolved_input / candidate)
            resolved_items.append(f"{key}={resolved_value}")
        parts.extend(f"override_input:{entry}" for entry in resolved_items)

    output_overrides = _parse_overrides(getattr(args, "override_output_stem", None))
    if output_overrides:
        parts.extend(
            f"override_output:{key}={value}"
            for key, value in sorted(output_overrides.items())
        )

    subcommand_overrides = _parse_overrides(
        getattr(args, "override_subcommand", None), allow_empty_value=True
    )
    if subcommand_overrides:
        parts.extend(
            f"override_subcommand:{key}={value}"
            for key, value in sorted(subcommand_overrides.items())
        )

    return "\n".join(parts)


def _hash_run_descriptor(descriptor: str) -> str:
    digest = hashlib.sha256(descriptor.encode("utf-8")).hexdigest()
    return digest[:32]


def _resolve_run_id(args: argparse.Namespace, *, descriptor: str) -> str:
    cli_value = getattr(args, "run_id", None)
    env_value = os.environ.get(_RUN_ID_ENV)
    for candidate in (cli_value, env_value):
        if candidate in (None, argparse.SUPPRESS):
            continue
        text = str(candidate).strip()
        if not text:
            raise ValueError("run identifier must not be empty")
        return text
    return _hash_run_descriptor(descriptor)


def _configure_logging(
    level_name: str,
    *,
    run_id: str | None = None,
    stream: IO[str] | None = None,
    handlers: Iterable[logging.Handler] | None = None,
) -> Logger:
    """Configure structured logging for the orchestration workflow."""

    normalised = level_name.upper()
    valid_levels = {"DEBUG", "INFO", "WARN", "WARNING", "ERROR"}
    if normalised not in valid_levels:
        raise ValueError(f"invalid log level: {level_name}")

    if run_id is None:
        raise ValueError("run_id must be provided")
    resolved_run_id = str(run_id).strip()
    if not resolved_run_id:
        raise ValueError("run_id must be provided")

    extra_handlers = list(handlers) if handlers is not None else []

    return configure_logger(
        LoggerConfig(
            level=normalised,
            run_id=resolved_run_id,
            stream=stream,
            handlers=extra_handlers,
        )
    )


def _parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Create the CLI parser and return parsed arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Run the document, target, assay, test item and activity pipelines "
            "sequentially with shared configuration."
        )
    )
    parser.add_argument(
        "--base-path",
        type=Path,
        default=Path("data"),
        help="Root directory containing input and output folders",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("input"),
        help="Directory with source CSV files relative to --base-path",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output"),
        help="Destination directory relative to --base-path",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="YAML configuration shared by all pipelines",
    )
    parser.add_argument(
        "--date",
        dest="date_prefix",
        default=None,
        help=(
            "Date prefix used to build output filenames. Defaults to "
            "local.io.default_date_prefix or CHEMBL_DA_DEFAULT_DATE_PREFIX"
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging verbosity for the orchestrator and child pipelines",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable debug logging for the orchestrator and delegated pipelines",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=(
            "Maximum number of identifiers processed by each pipeline; use 0 to "
            "skip execution"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run pipelines even if the target file already exists",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip pipeline execution when the output file is present",
    )
    parser.add_argument(
        "--rerun-postprocess",
        action="store_true",
        help=(
            "Rebuild stage-aligned exports even if previous runs already produced "
            "them"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Run pipelines without writing outputs or performing side effects, "
            "including output directory creation"
        ),
    )
    parser.add_argument(
        "--print-config",
        action="store_true",
        help="Print the resolved configuration and exit without running pipelines",
    )
    parser.add_argument(
        "--pipeline-registry",
        type=Path,
        default=None,
        help="Optional YAML file describing pipeline step definitions",
    )
    parser.add_argument(
        "--override-input",
        action="append",
        metavar="STEP=FILENAME",
        default=[],
        help="Override the input CSV filename for a pipeline step",
    )
    parser.add_argument(
        "--override-output-stem",
        action="append",
        metavar="STEP=STEM",
        default=[],
        help="Override the output filename stem for a pipeline step",
    )
    parser.add_argument(
        "--override-subcommand",
        action="append",
        metavar="STEP=SUBCOMMAND",
        default=[],
        help="Override the CLI subcommand used to invoke a pipeline step",
    )
    parser.add_argument(
        "--run-id",
        default=None,
        help="Override the computed run identifier",
    )

    parsed = parser.parse_args(argv)
    if not hasattr(parsed, "invocation"):
        parsed.invocation = resolve_invocation(parser.prog, argv)
    return parsed


def _prepare_config(
    args: argparse.Namespace, steps: Sequence[PipelineStep] | None = None
) -> PipelineRunConfig:
    """Validate CLI inputs and construct :class:`PipelineRunConfig`."""

    effective_steps = tuple(DEFAULT_PIPELINE_STEPS if steps is None else steps)
    input_files = {step.name: step.input_filename for step in effective_steps}
    output_stems = {step.name: step.output_stem for step in effective_steps}
    subcommands = {step.name: step.subcommand for step in effective_steps}

    base_path = args.base_path.expanduser().resolve()
    input_dir = _resolve_path(base_path, args.input_dir)
    output_dir = _resolve_path(base_path, args.output_dir)
    config_candidate = args.config or DEFAULT_CONFIG_PATH
    config_path = Path(config_candidate).expanduser().resolve()
    dry_run = bool(getattr(args, "dry_run", False))

    log_level = str(getattr(args, "log_level", "INFO")).upper()
    if getattr(args, "verbose", False):
        log_level = "DEBUG"

    if args.limit is not None and args.limit < 0:
        raise ValueError("--limit must be zero or a positive integer")

    if not input_dir.exists():
        raise FileNotFoundError(f"input directory not found: {input_dir}")
    if not config_path.exists():
        raise FileNotFoundError(f"configuration file not found: {config_path}")
    return PipelineRunConfig(
        base_path=base_path,
        input_dir=input_dir,
        output_dir=output_dir,
        config_path=config_path,
        date_prefix=args.date_prefix,
        log_level=log_level,
        limit=args.limit,
        force=args.force,
        skip_existing=args.skip_existing,
        dry_run=dry_run,
        rerun_postprocess=bool(getattr(args, "rerun_postprocess", False)),
        input_files=input_files,
        output_stems=output_stems,
        subcommands=subcommands,
    )


def _failure_sentinel_path(output_path: Path) -> Path:
    """Return the sentinel path recorded when a pipeline step fails."""

    return output_path.with_name(f"{output_path.name}.failed")


def _coerce_exit_code(value: object) -> int:
    """Translate ``SystemExit.code`` values into integer exit codes."""

    if isinstance(value, int):
        return value
    if value is None:
        return 0
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return 1
    return 1


@dataclass
class SidecarArtefact:
    """Describe auxiliary files associated with a pipeline output."""

    destination: Path
    final_path: Path | None = None
    working_path: Path | None = None


def _discover_sidecars(
    final_output: Path,
    working_output: Path,
    *,
    max_depth: int | None = None,
    include_patterns: Sequence[str] | None = None,
) -> dict[Path, SidecarArtefact]:
    """Return all auxiliary files derived from ``final_output`` and ``working_output``."""

    final_dir = final_output.parent
    working_dir = working_output.parent
    sentinel_name = _failure_sentinel_path(final_output).name
    patterns = {
        final_output.name,
        final_output.with_suffix("").name,
        working_output.name,
        working_output.with_suffix("").name,
    }
    prefix_candidates: set[str] = set()
    for candidate in (final_output, working_output):
        current = candidate
        while True:
            prefix_candidates.add(current.name)
            if not current.suffix:
                break
            current = current.with_suffix("")
    replacements = (
        (working_output.name, final_output.name),
        (working_output.with_suffix("").name, final_output.with_suffix("").name),
    )
    main_outputs = {final_output.resolve(), working_output.resolve()}

    def _matches_explicit_patterns(rel_path: Path, name: str) -> bool:
        if include_patterns is None:
            return False
        rel_value = rel_path.as_posix()
        return any(
            fnmatch(rel_value, pattern) or fnmatch(name, pattern)
            for pattern in include_patterns
        )

    def _normalise_relative(path: Path) -> Path:
        parts: list[str] = []
        for part in path.parts:
            normalised = part
            for old, new in replacements:
                if old == new or old not in normalised:
                    continue
                normalised = normalised.replace(old, new)
            parts.append(normalised)
        return Path(*parts)

    def _collect(base_dir: Path) -> dict[Path, Path]:
        collected: dict[Path, Path] = {}
        if not base_dir.exists():
            return collected
        queue: deque[tuple[Path, int, bool]] = deque([(base_dir, 0, False)])
        while queue:
            current, depth, within_branch = queue.popleft()
            try:
                entries = list(current.iterdir())
            except OSError:  # pragma: no cover - filesystem race protection
                continue
            for entry in entries:
                rel_path = entry.relative_to(base_dir)
                name = entry.name
                if entry.is_dir():
                    next_depth = depth + 1
                    if max_depth is not None and next_depth > max_depth:
                        continue
                    matches_prefix = any(
                        name.startswith(prefix) for prefix in prefix_candidates
                    )
                    matches_pattern = _matches_explicit_patterns(rel_path, name)
                    next_within = within_branch or matches_prefix or matches_pattern
                    if next_within:
                        queue.append((entry, next_depth, next_within))
                    continue
                matches_default = any(name.startswith(pattern) for pattern in patterns)
                matches_explicit = _matches_explicit_patterns(rel_path, name)
                if not matches_default and not matches_explicit:
                    continue
                if name == sentinel_name:
                    continue
                try:
                    resolved = entry.resolve()
                except OSError:  # pragma: no cover - filesystem race protection
                    continue
                if resolved in main_outputs:
                    continue
                collected[rel_path] = entry
        return collected

    sidecars: dict[Path, SidecarArtefact] = {}

    for rel_path, path in _collect(final_dir).items():
        canonical_rel = _normalise_relative(rel_path)
        destination = final_dir / canonical_rel
        entry = sidecars.get(canonical_rel)
        if entry is None:
            entry = SidecarArtefact(destination=destination)
            sidecars[canonical_rel] = entry
        entry.final_path = path

    for rel_path, path in _collect(working_dir).items():
        canonical_rel = _normalise_relative(rel_path)
        destination = final_dir / canonical_rel
        entry = sidecars.get(canonical_rel)
        if entry is None:
            entry = SidecarArtefact(destination=destination)
            sidecars[canonical_rel] = entry
        entry.working_path = path

    return sidecars


def _remove_path(path: Path) -> None:
    """Remove ``path`` handling transient Windows sharing violations."""

    if not path.exists():
        return
    for attempt in range(1, _UNLINK_MAX_ATTEMPTS + 1):
        try:
            path.unlink()
            return
        except FileNotFoundError:
            return
        except PermissionError as exc:
            if attempt == _UNLINK_MAX_ATTEMPTS:
                raise
            _LOGGER.debug(
                "unlink_retry_permission",
                path=str(path),
                attempt=attempt,
                error=str(exc),
                exc_info=exc,
            )
            time.sleep(_UNLINK_RETRY_SLEEP_SECONDS)
            continue
        except OSError as exc:
            winerror = getattr(exc, "winerror", None)
            if (
                winerror == _WINDOWS_SHARING_VIOLATION
                and attempt < _UNLINK_MAX_ATTEMPTS
            ):
                _LOGGER.debug(
                    "unlink_retry_sharing_violation",
                    path=str(path),
                    attempt=attempt,
                    error=str(exc),
                    exc_info=exc,
                )
                time.sleep(_UNLINK_RETRY_SLEEP_SECONDS)
                continue
            raise


def _cleanup_empty_directories(path: Path, *, root: Path) -> None:
    """Remove empty directories upward from ``path`` until reaching ``root``."""

    try:
        root_resolved = root.resolve()
    except OSError:  # pragma: no cover - path vanished concurrently
        return
    current = path
    while True:
        try:
            current_resolved = current.resolve()
        except FileNotFoundError:
            parent = current.parent
            if parent == current:
                break
            current = parent
            continue
        except OSError:  # pragma: no cover - defensive guard
            break
        if current_resolved == root_resolved or current == root:
            break
        try:
            current.rmdir()
        except OSError:
            break
        parent = current.parent
        if parent == current:
            break
        current = parent


def _run_step(
    step: PipelineStep,
    cfg: PipelineRunConfig,
    base_config: Config,
    input_path: Path,
    final_output: Path,
    working_output: Path,
) -> StepExecutionResult:
    """Execute ``step`` with ``cfg`` returning the resulting exit code."""

    if not input_path.exists():
        _LOGGER.error("step_input_missing", step=step.name, path=str(input_path))
        return StepExecutionResult(
            exit_code=1,
            executed=False,
            status="failed",
            reason="input_missing",
        )
    if cfg.skip_existing and final_output.exists() and not cfg.force:
        _LOGGER.info("step_skipped_existing", step=step.name, path=str(final_output))
        return StepExecutionResult(
            exit_code=0,
            executed=False,
            status="skipped",
            reason="skip_existing",
        )

    if cfg.limit == 0:
        _LOGGER.info("step_skip_limit", step=step.name, limit=cfg.limit)
        return StepExecutionResult(
            exit_code=0,
            executed=False,
            status="skipped",
            reason="limit",
        )

    api = _PIPELINE_APIS.get(step.name)
    if api is None:
        arguments = step.build_arguments(cfg, output_path=working_output)
        _LOGGER.debug("step_arguments", step=step.name, arguments=arguments)
        try:
            exit_code = step.main(arguments)
        except SystemExit as exc:
            exit_code = _coerce_exit_code(exc.code)
            status = "success" if exit_code == 0 else "failed"
            return StepExecutionResult(
                exit_code=exit_code,
                executed=True,
                status=status,
                reason="system_exit",
            )
        except BaseException:
            raise
        status = "success" if exit_code == 0 else "failed"
        return StepExecutionResult(
            exit_code=exit_code,
            executed=True,
            status=status,
            reason=None if status == "success" else "non_zero_exit",
        )

    options = api.build_options(cfg, input_path, working_output)
    result = api.runner(base_config, options)
    executed = bool(result.executed)
    if not executed and result.exit_code == 0:
        status = "skipped"
    else:
        status = "success" if result.exit_code == 0 else "failed"
    reason = result.reason
    if reason is None and status == "failed":
        reason = "non_zero_exit"
    return StepExecutionResult(
        exit_code=result.exit_code,
        executed=executed,
        status=status,
        reason=reason,
    )


def _finalize_step_success(
    final_output: Path, working_output: Path, sentinel_path: Path
) -> None:
    """Rename temporary outputs into place and clear failure sentinels."""

    sidecars = _discover_sidecars(final_output, working_output)
    working_dir = working_output.parent
    final_dir = final_output.parent

    if working_output.exists():
        if final_output.exists():
            _remove_path(final_output)
        working_output.replace(final_output)

    for sidecar in sidecars.values():
        working_path = sidecar.working_path
        if working_path is None or not working_path.exists():
            continue
        destination = sidecar.destination
        destination.parent.mkdir(parents=True, exist_ok=True)
        try:
            same_location = working_path.resolve() == destination.resolve()
        except OSError:
            same_location = False
        if same_location:
            final_path = sidecar.final_path
            if (
                final_path is not None
                and final_path.exists()
                and final_path != destination
            ):
                final_parent = final_path.parent
                _remove_path(final_path)
                _cleanup_empty_directories(final_parent, root=final_dir)
            continue
        if destination.exists():
            _remove_path(destination)
        original_parent = working_path.parent
        working_path.replace(destination)
        _cleanup_empty_directories(original_parent, root=working_dir)
        final_path = sidecar.final_path
        if final_path is not None and final_path.exists() and final_path != destination:
            final_parent = final_path.parent
            _remove_path(final_path)
            _cleanup_empty_directories(final_parent, root=final_dir)
    if sentinel_path.exists():
        _remove_path(sentinel_path)


def _resolve_postprocess_table(step: PipelineStep, final_output: Path) -> str | None:
    """Return the postprocess table identifier for ``step`` if supported."""

    table = getattr(step, "output_stem", "")
    if not table or table not in POSTPROCESS_SUPPORTED_TABLES:
        return None
    if final_output.suffix.lower() != ".csv":
        return None
    return str(table)


def _run_postprocess_hook(
    step: PipelineStep,
    final_output: Path,
    *,
    table: str | None = None,
) -> PostprocessResult | None:
    """Execute the post-processing pipeline for ``step`` when available."""

    if table is None:
        table = _resolve_postprocess_table(step, final_output)
    if table is None:
        return None
    if not final_output.exists():
        _LOGGER.warning(
            "postprocess_input_missing",
            step=step.name,
            table=table,
            input=str(final_output),
        )
        return None

    handlers = _POSTPROCESS_HANDLERS.get(table)
    if handlers is None:
        _LOGGER.warning(
            "postprocess_handlers_missing",
            step=step.name,
            table=table,
        )
        return None
    destination = final_output.with_name(f"output_postprocessed.{table}.csv")
    _LOGGER.info(
        "postprocess_start",
        step=step.name,
        table=table,
        input=str(final_output),
        output=str(destination),
    )
    pipeline_cfg = load_postprocess_pipeline_config(table, None)
    csv_cfg = get_postprocess_csv_config(pipeline_cfg)
    runtime_config = PostprocessingPipelineConfig(
        pipeline_config=pipeline_cfg,
        csv_runtime_config=csv_cfg,
        runner=handlers.runner,
        validator=handlers.validator,
        schema=handlers.schema,
        logger=_LOGGER,
    )
    pipeline_result = run_postprocessing_pipeline(
        table,
        final_output,
        destination,
        runtime_config,
    )
    metrics = pipeline_result.metrics
    report_path = pipeline_result.report_path
    if metrics is None:
        raise RuntimeError(f"postprocess metrics missing for table {table}")
    if report_path is None:
        raise RuntimeError(f"postprocess report missing for table {table}")
    _LOGGER.info(
        "postprocess_done",
        step=step.name,
        table=table,
        output=str(pipeline_result.output_path),
        report=str(report_path),
        rows=int(metrics.output_rows),
        columns=int(metrics.output_columns),
    )
    return PostprocessResult(
        table=table,
        output_path=pipeline_result.output_path,
        report_path=report_path,
        metrics=metrics,
    )


def _cleanup_failed_step(
    final_output: Path,
    working_output: Path,
    sentinel_path: Path,
    *,
    executed: bool,
) -> None:
    """Remove partial outputs and persist a failure sentinel."""

    sidecars = _discover_sidecars(final_output, working_output)
    working_dir = working_output.parent
    final_dir = final_output.parent

    candidates = [working_output]
    if executed:
        candidates.append(final_output)
    for candidate in candidates:
        if candidate.exists():
            try:
                _remove_path(candidate)
            except OSError as exc:  # pragma: no cover - defensive guard
                _LOGGER.warning(
                    "step_cleanup_failed",
                    path=str(candidate),
                    error=str(exc),
                )

    for sidecar in sidecars.values():
        destination = sidecar.destination
        is_failure = destination.name.endswith("_failure_cases.csv")
        working_path = sidecar.working_path
        final_path = sidecar.final_path

        if is_failure:
            if working_path is not None and working_path.exists():
                destination.parent.mkdir(parents=True, exist_ok=True)
                try:
                    if destination.exists():
                        _remove_path(destination)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(destination),
                        error=str(exc),
                    )
                try:
                    original_parent = working_path.parent
                    working_path.replace(destination)
                    _cleanup_empty_directories(original_parent, root=working_dir)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_promote_failed",
                        path=str(working_path),
                        error=str(exc),
                    )
            elif executed and destination.exists():
                try:
                    _remove_path(destination)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(destination),
                        error=str(exc),
                    )
                else:
                    _cleanup_empty_directories(destination.parent, root=final_dir)
            if (
                final_path is not None
                and final_path.exists()
                and final_path != destination
            ):
                try:
                    final_parent = final_path.parent
                    _remove_path(final_path)
                    _cleanup_empty_directories(final_parent, root=final_dir)
                except OSError as exc:  # pragma: no cover - defensive guard
                    _LOGGER.warning(
                        "failure_sidecar_cleanup_failed",
                        path=str(final_path),
                        error=str(exc),
                    )
            continue

        if working_path is not None and working_path.exists():
            try:
                original_parent = working_path.parent
                _remove_path(working_path)
                _cleanup_empty_directories(original_parent, root=working_dir)
            except OSError as exc:  # pragma: no cover - defensive guard
                _LOGGER.warning(
                    "step_cleanup_failed",
                    path=str(working_path),
                    error=str(exc),
                )
        if executed:
            removal_targets = {destination}
            if final_path is not None:
                removal_targets.add(final_path)
            for path in removal_targets:
                if path.exists():
                    try:
                        parent = path.parent
                        _remove_path(path)
                        _cleanup_empty_directories(parent, root=final_dir)
                    except OSError as exc:  # pragma: no cover - defensive guard
                        _LOGGER.warning(
                            "step_cleanup_failed",
                            path=str(path),
                            error=str(exc),
                        )
    try:
        sentinel_path.touch()
    except OSError as exc:  # pragma: no cover - defensive guard
        _LOGGER.warning(
            "sentinel_write_failed",
            path=str(sentinel_path),
            error=str(exc),
        )


def _warm_parent_catalog(cfg: PipelineRunConfig, base_config: Config) -> None:
    """Ensure the molecule parent catalogue cache exists before test item runs."""

    from library.integration.molecule_catalog import load_parent_catalog

    start_time = time.perf_counter()
    chembl_sources = base_config.sources.chembl
    catalog_cfg = chembl_sources.molecule_catalog
    cache_path = catalog_cfg.cache_path
    sqlite_path = catalog_cfg.sqlite_path

    cache_ready = cache_path.is_file()
    sqlite_ready = sqlite_path.is_file()

    log_kwargs: dict[str, str] = {
        "cache": str(cache_path),
        "sqlite": str(sqlite_path),
    }

    if cache_ready and sqlite_ready:
        elapsed = time.perf_counter() - start_time
        _LOGGER.info("parent_catalog_warm_skip", elapsed=elapsed, **log_kwargs)
        return

    testitem_cfg = chembl_sources.pipelines.testitem
    _LOGGER.info("parent_catalog_warm_start", **log_kwargs)

    def _catalog_client_factory(context: ETLContext) -> ChemblClient:
        return ChemblClient(
            api=chembl_sources.api,
            retry=base_config.system.retry,
            chembl=chembl_sources.cache,
            global_limiter=context.global_limiter,
        )

    try:
        with ETLContext(
            base_config, chembl_client_factory=_catalog_client_factory
        ) as context:
            load_parent_catalog(
                client=context.chembl_client,
                api_cfg=chembl_sources.api,
                catalog_cfg=catalog_cfg,
                timeout=testitem_cfg.timeout,
            )
    except TimeoutError as exc:
        elapsed = time.perf_counter() - start_time
        _LOGGER.error(
            "parent_catalog_warm_failed",
            elapsed=elapsed,
            error=str(exc),
            exc_info=exc,
            **log_kwargs,
        )
        raise
    except Exception as exc:  # pragma: no cover - defensive guard
        elapsed = time.perf_counter() - start_time
        context: dict[str, object] = {
            **log_kwargs,
            "elapsed": elapsed,
            "error": str(exc),
        }
        _LOGGER.exception(
            "parent_catalog_warm_failed",
            exc=exc,
            **context,
        )
        raise
    elapsed = time.perf_counter() - start_time
    _LOGGER.info("parent_catalog_warm_done", elapsed=elapsed, **log_kwargs)


def _compute_file_checksum(path: Path, *, chunk_size: int = 65536) -> str:
    """Return the hexadecimal SHA256 checksum for ``path``."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _describe_file(path: Path) -> dict[str, Any]:
    """Return manifest metadata for ``path`` including checksum when available."""

    info: dict[str, Any] = {"path": str(path)}
    try:
        exists = path.exists()
    except OSError:
        exists = False
    if not exists or not path.is_file():
        info.update(
            {
                "exists": False,
                "size_bytes": None,
                "checksum_sha256": None,
            }
        )
        return info
    try:
        stat_result = path.stat()
    except OSError:
        info.update(
            {
                "exists": False,
                "size_bytes": None,
                "checksum_sha256": None,
            }
        )
        return info
    info.update(
        {
            "exists": True,
            "size_bytes": stat_result.st_size,
            "checksum_sha256": _compute_file_checksum(path),
        }
    )
    return info


def _describe_sidecars(
    final_output: Path, working_output: Path
) -> list[dict[str, Any]]:
    """Return manifest metadata for sidecars associated with ``final_output``."""

    include_patterns = [
        f"*{final_output.name}",
        f"*{final_output.with_suffix('').name}",
        f"*{working_output.name}",
        f"*{working_output.with_suffix('').name}",
    ]
    stem = final_output.stem
    table_candidate: str | None = None
    if stem.startswith("output."):
        remainder = stem[len("output.") :]
        if "_" in remainder:
            table_candidate = remainder.split("_", 1)[0]
        elif remainder:
            table_candidate = remainder
    if table_candidate:
        include_patterns.extend(
            [
                f"*output_postprocessed.{table_candidate}.csv",
                f"*{table_candidate}.postprocess.report.json",
            ]
        )
    sidecars = _discover_sidecars(
        final_output,
        working_output,
        include_patterns=tuple(include_patterns),
    )
    described: list[dict[str, Any]] = []
    for artefact in sorted(sidecars.values(), key=lambda item: str(item.destination)):
        candidates: tuple[Path | None, ...] = (
            artefact.destination,
            artefact.final_path,
            artefact.working_path,
        )
        selected: Path | None = None
        for candidate in candidates:
            if candidate is None:
                continue
            if candidate.exists():
                selected = candidate
                break
        if selected is None:
            selected = artefact.destination
        described.append(_describe_file(selected))
    return described


def _complete_manifest_entry(
    entry: dict[str, Any],
    *,
    final_output: Path,
    working_output: Path,
    started_at: float,
) -> None:
    """Populate completion metadata for a manifest ``entry``."""

    entry["completed_at"] = datetime.now(UTC).isoformat()
    entry["duration_sec"] = round(time.perf_counter() - started_at, 6)
    entry["output"] = _describe_file(final_output)
    entry["sidecars"] = _describe_sidecars(final_output, working_output)
    report = load_output_report(final_output)
    if report is not None:
        merge_run_output(entry, report)


def _pending_manifest_entry(
    step: PipelineStep, cfg: PipelineRunConfig
) -> dict[str, Any]:
    """Return a manifest entry for ``step`` that has not been executed."""

    final_output = step.expected_output(cfg)
    working_output = temporary_output_path(final_output)
    entry: dict[str, Any] = {
        "name": step.name,
        "status": "pending",
        "exit_code": None,
        "executed": False,
        "reason": None,
        "started_at": None,
        "completed_at": None,
        "duration_sec": None,
        "output": _describe_file(final_output),
        "sidecars": _describe_sidecars(final_output, working_output),
    }
    return entry


def _write_run_manifest(
    cfg: PipelineRunConfig,
    *,
    run_started_at: datetime,
    run_completed_at: datetime,
    duration_seconds: float,
    exit_code: int,
    steps: Sequence[dict[str, Any]],
    run_id: str | None = None,
) -> None:
    """Persist the manifest for the pipeline execution."""

    manifest = {
        "run": {
            "started_at": run_started_at.isoformat(),
            "completed_at": run_completed_at.isoformat(),
            "duration_sec": round(duration_seconds, 6),
            "exit_code": exit_code,
            "status": "success" if exit_code == 0 else "failed",
            "date_prefix": cfg.date_prefix,
            "base_path": str(cfg.base_path),
            "input_dir": str(cfg.input_dir),
            "output_dir": str(cfg.output_dir),
            "config_path": str(cfg.config_path),
            "log_level": cfg.log_level,
            "limit": cfg.limit,
            "force": cfg.force,
            "skip_existing": cfg.skip_existing,
            "dry_run": cfg.dry_run,
        },
        "steps": list(steps),
    }

    if run_id is not None:
        manifest["run"]["run_id"] = run_id

    reports_dir = cfg.base_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    timestamp = run_started_at.astimezone(UTC).strftime("%Y%m%dT%H%M%S%fZ")
    manifest_path = reports_dir / f"run_{timestamp}.json"
    try:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
    except OSError as exc:  # pragma: no cover - defensive guard
        _LOGGER.warning(
            "manifest_write_failed",
            path=str(manifest_path),
            error=str(exc),
        )
        return

    latest_alias = reports_dir / "run_manifest.json"
    alias_removed = True
    if latest_alias.exists() or latest_alias.is_symlink():
        try:
            latest_alias.unlink()
        except OSError as exc:  # pragma: no cover - defensive guard
            alias_removed = False
            _LOGGER.warning(
                "manifest_alias_cleanup_failed",
                path=str(latest_alias),
                error=str(exc),
            )

    if not alias_removed:
        alias_tmp = latest_alias.with_name(f"{latest_alias.name}.tmp")
        try:
            alias_tmp.write_bytes(manifest_path.read_bytes())
            os.replace(alias_tmp, latest_alias)
        except OSError as exc:
            try:
                alias_tmp.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort cleanup
                pass
            _LOGGER.error(
                "manifest_alias_fallback_failed",
                path=str(latest_alias),
                temp_path=str(alias_tmp),
                target=str(manifest_path),
                error=str(exc),
            )
            raise RuntimeError("failed to update run manifest alias") from exc
        return

    try:
        latest_alias.symlink_to(manifest_path.name)
    except (OSError, NotImplementedError):  # pragma: no cover - platform dependent
        try:
            shutil.copy2(manifest_path, latest_alias)
        except OSError as exc:  # pragma: no cover - defensive guard
            _LOGGER.warning(
                "manifest_alias_update_failed",
                path=str(latest_alias),
                target=str(manifest_path),
                error=str(exc),
            )


@dataclass(frozen=True)
class PipelineExecutionPlan:
    """Describe the resolved execution order and artefact dependencies."""

    steps: tuple[PipelineStep, ...]
    dependencies: Mapping[str, frozenset[str]]
    produced_by: Mapping[str, str]
    external_artifacts: Mapping[str, tuple[str, ...]]


def _build_execution_plan(
    steps: Sequence[PipelineStep],
) -> PipelineExecutionPlan:
    """Return a deterministic execution plan for ``steps``."""

    if not steps:
        empty_mapping: dict[str, tuple[str, ...]] = {}
        return PipelineExecutionPlan(
            steps=(),
            dependencies={},
            produced_by={},
            external_artifacts=empty_mapping,
        )

    by_name: dict[str, PipelineStep] = {}
    for step in steps:
        if step.name in by_name:
            raise ValueError(f"duplicate pipeline step name: {step.name}")
        by_name[step.name] = step

    produced_by: dict[str, str] = {}
    for step in steps:
        produced = step.produces or (step.output_stem,)
        for artefact in produced:
            current = produced_by.get(artefact)
            if current is not None and current != step.name:
                raise ValueError(f"artefact '{artefact}' declared by multiple steps")
            produced_by[artefact] = step.name

    dependencies: dict[str, set[str]] = {
        step.name: set(step.depends_on) for step in steps
    }
    for step in steps:
        for artefact in step.consumes:
            producer = produced_by.get(artefact)
            if producer is not None and producer != step.name:
                dependencies[step.name].add(producer)

    missing: list[str] = []
    for name, deps in dependencies.items():
        unknown = sorted(dep for dep in deps if dep not in by_name)
        if unknown:
            missing.append(f"{name}: {', '.join(unknown)}")
    if missing:
        raise ValueError(
            "pipeline references unknown dependency: " + "; ".join(missing)
        )

    adjacency: dict[str, set[str]] = {name: set() for name in by_name}
    indegree: dict[str, int] = {}
    for name, deps in dependencies.items():
        indegree[name] = len(deps)
        for dep in deps:
            adjacency[dep].add(name)

    order_index = {step.name: index for index, step in enumerate(steps)}
    queue: list[tuple[int, str]] = []
    for name, degree in indegree.items():
        if degree == 0:
            heappush(queue, (order_index[name], name))

    ordered: list[str] = []
    while queue:
        _, current = heappop(queue)
        ordered.append(current)
        neighbours = sorted(adjacency[current], key=order_index.__getitem__)
        for neighbour in neighbours:
            indegree[neighbour] -= 1
            if indegree[neighbour] == 0:
                heappush(queue, (order_index[neighbour], neighbour))

    if len(ordered) != len(steps):
        remaining = sorted(name for name, degree in indegree.items() if degree > 0)
        raise ValueError("cyclic pipeline dependency detected: " + ", ".join(remaining))

    ordered_steps = tuple(by_name[name] for name in ordered)
    scheduled = set(ordered)
    external: dict[str, tuple[str, ...]] = {}
    for step in ordered_steps:
        requirements: list[str] = []
        for artefact in step.consumes:
            producer = produced_by.get(artefact)
            if producer is None or producer not in scheduled:
                requirements.append(artefact)
        external[step.name] = tuple(requirements)

    frozen_dependencies = {name: frozenset(deps) for name, deps in dependencies.items()}
    return PipelineExecutionPlan(
        steps=ordered_steps,
        dependencies=frozen_dependencies,
        produced_by=dict(produced_by),
        external_artifacts=external,
    )


def run_pipeline(
    cfg: PipelineRunConfig,
    *,
    steps: Sequence[PipelineStep] | None = None,
) -> int:
    """Execute all configured steps and return the resulting exit status."""

    effective_steps = tuple(DEFAULT_PIPELINE_STEPS if steps is None else steps)
    try:
        plan = _build_execution_plan(effective_steps)
    except ValueError as exc:
        _LOGGER.error("pipeline_schedule_error", error=str(exc), exc_info=exc)
        return 1
    effective_steps = plan.steps
    external_requirements = plan.external_artifacts
    try:
        base_config = load_config(cfg.config_path, base_path=cfg.base_path)
    except (ConfigError, ConfigLoaderError, ValidationError) as exc:
        _LOGGER.error("config_load_failed", error=str(exc), exc_info=exc)
        _LOGGER.info("pipeline_done", stage="pipeline", exit_code=1)
        return 1
    try:
        if not cfg.dry_run:
            ensure_dirs(base_config)
    except (FileNotFoundError, NotADirectoryError, OSError) as exc:
        _LOGGER.error("config_directory_error", error=str(exc), exc_info=exc)
        _LOGGER.info("pipeline_done", stage="pipeline", exit_code=1)
        return 1
    overall_status = 0
    run_started_at = datetime.now(UTC)
    run_started_clock = time.perf_counter()
    manifest_entries: list[dict[str, Any]] = []
    failed_index: int | None = None
    last_executed_index = -1
    completed_steps: set[str] = set()
    _LOGGER.info("pipeline_start", stage="pipeline")

    def _prepare_step(
        step: PipelineStep,
        entry: dict[str, Any],
        index: int,
    ) -> PreparedPipelineStep:
        def _invoke(
            current_cfg: PipelineRunConfig,
            input_path: Path,
            final_output: Path,
            working_output: Path,
        ) -> StepExecutionResult:
            nonlocal overall_status, failed_index, last_executed_index

            sentinel_path = _failure_sentinel_path(final_output)
            entry["started_at"] = datetime.now(UTC).isoformat()
            step_started_clock = time.perf_counter()
            _LOGGER.info("step_start", step=step.name)

            if current_cfg.dry_run:
                _LOGGER.info("step_skip_dry_run", step=step.name)
                entry.update(
                    {
                        "status": "skipped",
                        "exit_code": 0,
                        "executed": False,
                        "reason": "dry_run",
                    }
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                last_executed_index = index
                completed_steps.add(step.name)
                return StepExecutionResult(
                    exit_code=0,
                    executed=False,
                    status="skipped",
                    reason="dry_run",
                )

            missing_external: list[tuple[str, Path]] = []
            required_external = external_requirements.get(step.name, ())
            if required_external:
                for artefact in required_external:
                    candidate = _resolve_consumed_artifact_path(current_cfg, artefact)
                    if not candidate.exists():
                        missing_external.append((artefact, candidate))
            if missing_external:
                overall_status = 1
                entry.update(
                    {
                        "status": "failed",
                        "exit_code": overall_status,
                        "executed": False,
                        "reason": "dependency_missing",
                    }
                )
                _LOGGER.error(
                    "step_dependencies_missing",
                    step=step.name,
                    missing=[
                        {"artefact": artefact, "path": str(path)}
                        for artefact, path in missing_external
                    ],
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                failed_index = index
                last_executed_index = index
                return StepExecutionResult(
                    exit_code=1,
                    executed=False,
                    status="failed",
                    reason="dependency_missing",
                )

            if working_output.exists():
                _remove_path(working_output)

            if step.name == "testitem":
                should_skip_warm = (
                    current_cfg.skip_existing
                    and final_output.exists()
                    and not current_cfg.force
                )
                if not should_skip_warm:
                    try:
                        _warm_parent_catalog(current_cfg, base_config)
                    except TimeoutError as exc:
                        overall_status = 1
                        entry.update(
                            {
                                "status": "failed",
                                "exit_code": overall_status,
                                "executed": False,
                                "reason": "parent_catalog_timeout",
                            }
                        )
                        _LOGGER.error(
                            "parent_catalog_warm_timeout", error=str(exc), exc_info=exc
                        )
                        _complete_manifest_entry(
                            entry,
                            final_output=final_output,
                            working_output=working_output,
                            started_at=step_started_clock,
                        )
                        failed_index = index
                        last_executed_index = index
                        return StepExecutionResult(
                            exit_code=1,
                            executed=False,
                            status="failed",
                            reason="parent_catalog_timeout",
                        )
                    except Exception as exc:  # pragma: no cover - defensive guard
                        overall_status = 1
                        entry.update(
                            {
                                "status": "failed",
                                "exit_code": overall_status,
                                "executed": False,
                                "reason": "parent_catalog_error",
                            }
                        )
                        _LOGGER.error(
                            "parent_catalog_warm_error", error=str(exc), exc_info=exc
                        )
                        _complete_manifest_entry(
                            entry,
                            final_output=final_output,
                            working_output=working_output,
                            started_at=step_started_clock,
                        )
                        failed_index = index
                        last_executed_index = index
                        return StepExecutionResult(
                            exit_code=1,
                            executed=False,
                            status="failed",
                            reason="parent_catalog_error",
                        )

            try:
                result = _run_step(
                    step,
                    current_cfg,
                    base_config,
                    input_path,
                    final_output,
                    working_output,
                )
            except SystemExit as exc:  # pragma: no cover - defensive guard
                exit_code = _coerce_exit_code(exc.code)
                entry.update(
                    {
                        "status": "failed",
                        "exit_code": exit_code,
                        "executed": True,
                        "reason": "system_exit",
                    }
                )
                _LOGGER.error("step_system_exit", step=step.name, exit_code=exit_code)
                _cleanup_failed_step(
                    final_output,
                    working_output,
                    sentinel_path,
                    executed=True,
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                overall_status = exit_code
                failed_index = index
                last_executed_index = index
                return StepExecutionResult(
                    exit_code=exit_code,
                    executed=True,
                    status="failed",
                    reason="system_exit",
                )
            except BaseException as exc:  # pragma: no cover - defensive guard
                entry.update(
                    {
                        "status": "failed",
                        "exit_code": 1,
                        "executed": True,
                        "reason": "exception",
                    }
                )
                _LOGGER.exception("step_exception", step=step.name, error=str(exc))
                _cleanup_failed_step(
                    final_output,
                    working_output,
                    sentinel_path,
                    executed=True,
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                overall_status = 1
                failed_index = index
                last_executed_index = index
                return StepExecutionResult(
                    exit_code=1,
                    executed=True,
                    status="failed",
                    reason="exception",
                )

            entry.update(
                {
                    "exit_code": result.exit_code,
                    "executed": result.executed,
                    "reason": result.reason,
                }
            )

            last_executed_index = index

            if result.exit_code != 0:
                _LOGGER.error("step_failed", step=step.name, exit_code=result.exit_code)
                entry["status"] = result.status
                _cleanup_failed_step(
                    final_output,
                    working_output,
                    sentinel_path,
                    executed=result.executed,
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                overall_status = result.exit_code
                failed_index = index
                return result

            entry["status"] = result.status
            output_existed = final_output.exists() or working_output.exists()
            _finalize_step_success(final_output, working_output, sentinel_path)
            if (
                result.executed
                and result.exit_code == 0
                and not output_existed
                and not final_output.exists()
            ):
                entry.update(
                    {
                        "status": "failed",
                        "exit_code": 1,
                        "executed": True,
                        "reason": "output_missing",
                    }
                )
                _LOGGER.error(
                    "step_output_missing",
                    step=step.name,
                    path=str(final_output),
                )
                _cleanup_failed_step(
                    final_output,
                    working_output,
                    sentinel_path,
                    executed=True,
                )
                _complete_manifest_entry(
                    entry,
                    final_output=final_output,
                    working_output=working_output,
                    started_at=step_started_clock,
                )
                overall_status = 1
                failed_index = index
                last_executed_index = index
                return StepExecutionResult(
                    exit_code=1,
                    executed=True,
                    status="failed",
                    reason="output_missing",
                )
            postprocess_result: PostprocessResult | None = None
            postprocess_table = _resolve_postprocess_table(step, final_output)
            if result.executed and postprocess_table is not None:
                try:
                    postprocess_result = _run_postprocess_hook(
                        step,
                        final_output,
                        table=postprocess_table,
                    )
                except Exception as exc:  # pragma: no cover - defensive guard
                    entry.update(
                        {
                            "status": "failed",
                            "exit_code": 1,
                            "reason": "postprocess_failed",
                        }
                    )
                    _LOGGER.exception(
                        "postprocess_failed",
                        step=step.name,
                        table=postprocess_table,
                        input=str(final_output),
                        error=str(exc),
                    )
                    _cleanup_failed_step(
                        final_output,
                        working_output,
                        sentinel_path,
                        executed=True,
                    )
                    _complete_manifest_entry(
                        entry,
                        final_output=final_output,
                        working_output=working_output,
                        started_at=step_started_clock,
                    )
                    overall_status = 1
                    failed_index = index
                    return StepExecutionResult(
                        exit_code=1,
                        executed=True,
                        status="failed",
                        reason="postprocess_failed",
                    )
            if postprocess_result is not None:
                entry["postprocess"] = {
                    "output": str(postprocess_result.output_path),
                    "report": str(postprocess_result.report_path),
                    "rows": postprocess_result.metrics.output_rows,
                    "columns": postprocess_result.metrics.output_columns,
                }
            _complete_manifest_entry(
                entry,
                final_output=final_output,
                working_output=working_output,
                started_at=step_started_clock,
            )
            _LOGGER.info("step_done", step=step.name)
            completed_steps.add(step.name)
            last_executed_index = index
            return result

        return PreparedPipelineStep(step=step, invoke=_invoke)

    prepared_steps: list[PreparedPipelineStep] = []
    for index, step in enumerate(effective_steps):
        entry = _pending_manifest_entry(step, cfg)
        manifest_entries.append(entry)
        prepared_steps.append(_prepare_step(step, entry, index))

    for _ in execute_workflow(cfg, prepared_steps):
        pass

    if failed_index is not None:
        for pending_entry in manifest_entries[failed_index + 1 :]:
            pending_entry.update(
                {
                    "status": "blocked",
                    "executed": False,
                    "reason": "dependency_failed",
                }
            )

    if failed_index is None and last_executed_index + 1 == len(effective_steps):
        _LOGGER.info("workflow_complete")

    run_completed_at = datetime.now(UTC)
    duration_seconds = time.perf_counter() - run_started_clock
    logger_cfg = getattr(_LOGGER, "_cfg", None)
    manifest_run_id = getattr(logger_cfg, "run_id", None)
    _write_run_manifest(
        cfg,
        run_started_at=run_started_at,
        run_completed_at=run_completed_at,
        duration_seconds=duration_seconds,
        exit_code=overall_status,
        steps=manifest_entries,
        run_id=manifest_run_id,
    )

    _LOGGER.info("pipeline_done", stage="pipeline", exit_code=overall_status)
    return overall_status


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for the orchestration script."""

    args = _parse_args(argv)
    desired_level = (
        "DEBUG" if getattr(args, "verbose", False) else str(args.log_level).upper()
    )
    if desired_level not in {"DEBUG", "INFO", "WARN", "WARNING", "ERROR"}:
        raise SystemExit(f"invalid log level: {args.log_level}")

    script_name = Path(__file__).with_suffix("").name
    resolved_base_path = Path(args.base_path).expanduser().resolve()
    try:
        _ensure_date_prefix(args, base_path=resolved_base_path)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    try:
        descriptor = _canonical_run_descriptor(args, base_path=resolved_base_path)
        resolved_run_id = _resolve_run_id(args, descriptor=descriptor)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    base_cfg = LoggerConfig(level=desired_level, run_id=resolved_run_id)
    log_directory = resolved_base_path / "logs"
    status = 1

    with setup_cli_logging(
        script_name,
        base_cfg,
        getattr(args, "date_prefix", None),
        log_dir=log_directory,
    ) as logging_ctx:
        try:
            logger = _configure_logging(
                logging_ctx.log_cfg.level,
                run_id=logging_ctx.log_cfg.run_id,
                stream=logging_ctx.console_stream,
                handlers=logging_ctx.log_cfg.handlers,
            )
        except ValueError as exc:
            raise SystemExit(str(exc)) from exc

        global _LOGGER
        _LOGGER = logger

        try:
            steps = _resolve_pipeline_steps(args)
        except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
            _LOGGER.error("registry_error", error=str(exc), exc_info=exc)
            status = 1
        else:
            try:
                cfg = _prepare_config(args, steps)
            except (FileNotFoundError, OSError, ValueError) as exc:
                _LOGGER.error("configuration_error", error=str(exc), exc_info=exc)
                status = 1
            else:
                if getattr(args, "print_config", False):
                    try:
                        config_obj = load_config(
                            cfg.config_path, base_path=cfg.base_path
                        )
                    except (ConfigError, ConfigLoaderError, ValidationError) as exc:
                        _LOGGER.error(
                            "config_load_failed", error=str(exc), exc_info=exc
                        )
                        status = 1
                    else:
                        print_config(config_obj)
                        status = 0
                else:
                    status = run_pipeline(cfg, steps=steps)
                    if status != 0:
                        _LOGGER.error("workflow_failed", exit_code=status)
                    else:
                        _LOGGER.info("workflow_succeeded")

    return status


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())


def _temporary_output_path(output_path: Path) -> Path:
    """Backward compatible alias for :func:`temporary_output_path`."""

    return temporary_output_path(output_path)
