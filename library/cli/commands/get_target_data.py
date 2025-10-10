"""Command line interface for retrieving target data from external sources.

Example
-------
Fetch ChEMBL target information for identifiers in ``targets.csv``::

    chembl-data-acquisition get-target-data chembl --config config/config.yaml --input targets.csv
"""

# Changelog:
# - Batch column assignments in merge logic to avoid pandas fragmentation warnings.

from __future__ import annotations

import argparse
import errno
import math
import os
import re
import shutil
import stat
import sys
import time
from datetime import datetime, timezone
from collections.abc import Collection, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass, replace
from itertools import islice
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd
import requests
from pandera.errors import SchemaErrors

import library.cli_utils as cli_utils_module
from library import SidecarErrors, cli, io
from library.cli import (
    Logger,
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
)
from library.cli.base import PipelineCLIBase
from library.cli.logging import CLILoggingContext
from library.cli.pipeline_definition import normalise_definition
from library.cli_utils import (
    PipelineError,
    PipelineExecutionResult,
    run_cli_command,
    run_pipeline,
)
from library.common.csv_utils import write_csv_deterministic
from library.common.log import logger
from library.config import (
    Config,
    _serialize_paths,
)
from library.integration import chembl_library as cl
from library.integration import iuphar_library as ii
from library.integration import uniprot_library as uu
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.orchestration import ETLContext
from library.pipelines.common import PipelineRunResult, add_pipeline_metadata
from library.pipelines.common.metadata import get_pipeline_version
from library.pipelines.target import TargetPipelineOptions
from library.pipelines.target import postprocessing as tp
from library.pipelines.target import protein_classification as pc
from library.pipelines.target.chembl_target import normalize_reaction_ec_numbers
from library.pipelines.target.defaults import TARGET_MODE_DEFAULTS, ModeDefaults
from library.qa.reporting import build_table_quality_hook, is_quality_enabled
from library.utils.data_correlation import generate_correlation_report
from library.utils.qc_report import generate_qc_report
from library.schemas import TargetsSchema, normalize_targets
from library.schemas.targets import TARGETS_COLUMN_ORDER
from library.validation import ValidationResult, validate_targets

if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from library.postprocess import PostprocessingPipelineResult

@contextmanager
def _override_cli_meta_writer() -> Iterator[None]:
    """Temporarily patch CLI metadata writer used by ``run_pipeline``."""

    original_cli_write_meta = cli_utils_module.write_meta_yaml
    cli_utils_module.write_meta_yaml = write_meta_yaml
    try:
        yield
    finally:
        cli_utils_module.write_meta_yaml = original_cli_write_meta


def _run_pipeline_with_meta(**kwargs: object) -> PipelineExecutionResult:
    """Invoke :func:`run_pipeline` with project-specific metadata writer."""

    with _override_cli_meta_writer():
        params = dict(kwargs)
        try:
            fetcher = params.pop("fetcher")
            output_path = params.pop("output_path")
            failure_path = params.pop("failure_path")
        except KeyError as exc:  # pragma: no cover - defensive validation
            missing = exc.args[0]
            raise TypeError(
                f"run_pipeline missing required argument: {missing}"
            ) from exc

        cfg = params.pop("cfg", None)
        logger = params.pop("logger", None)
        definition = params.pop("definition", None)

        pipeline_definition = normalise_definition(definition, params)

        return run_pipeline(
            definition=pipeline_definition,
            fetcher=fetcher,
            output_path=output_path,
            failure_path=failure_path,
            cfg=cfg,
            logger=logger,
        )


UNIPROT_MISSING_VALUE = ""


DEFAULT_INPUT_NAME = "target.csv"
DEFAULT_OUTPUT_STEM = "targets"
RAW_SUFFIX = "_raw"
NORMALIZED_SUFFIX = "_normalized"
COMMAND_CHOICES: tuple[str, ...] = ("uniprot", "chembl", "iuphar", "all")


_IUPHAR_OVERRIDE_COLUMNS: frozenset[str] = frozenset(
    {
        "target_id",
        "GuidetoPHARMACOLOGY",
        "gtop_synonyms",
        "gtop_natural_ligands_n",
        "gtop_interactions_n",
        "gtop_function_text_short",
        "IUPHAR_family_id",
        "IUPHAR_type",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "IUPHAR_chain",
        "IUPHAR_name",
        "iuphar_name",
        "full_id_path",
        "full_name_path",
        "iuphar_synonyms",
    }
)


def _matches_expected_target_input_name(filename: str) -> bool:
    """Return ``True`` if *filename* matches canonical target export patterns."""

    try:
        from library.postprocessing.target import _matches_expected_input_name
    except ImportError:  # pragma: no cover - optional dependency
        return False

    try:
        return bool(_matches_expected_input_name(filename))
    except Exception:  # pragma: no cover - defensive guard
        return False


def _prepare_targets_for_schema(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str], set[str]]:
    """Align *frame* to the :class:`TargetsSchema` structure lazily."""

    from library.postprocessing.target.export import (  # type: ignore
        prepare_targets_for_schema as _prepare,
    )

    return _prepare(frame)


class StoreWithSource(argparse.Action):
    """Store CLI values while tracking their origin."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: object,
        option_string: str | None = None,
    ) -> None:
        overrides = getattr(namespace, "_cli_overrides", None)
        if overrides is None:
            overrides = set()
            namespace._cli_overrides = overrides
        overrides.add(self.dest)
        setattr(namespace, self.dest, values)


def _non_negative_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:  # pragma: no cover - delegated to argparse
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or positive")
    return parsed


def _positive_float(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:  # pragma: no cover - delegated to argparse
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return parsed


@dataclass
class ParameterLogEntry:
    scope: str
    name: str
    value: object
    source: str


@dataclass
class IsoformPostprocessContext:
    args: argparse.Namespace | None = None
    http_requests: int | None = None


def _bool_from_cli(value: object) -> bool | None:
    """Return boolean representation for CLI flags, tolerating ``None``."""

    if value in (None, argparse.SUPPRESS):
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"", "none", "null"}:
            return None
        if lowered in {"1", "true", "yes", "on"}:
            return True
        if lowered in {"0", "false", "no", "off"}:
            return False
    return bool(value)


_TEMPORARY_EXPORT_SUFFIXES: tuple[str, ...] = (".tmp",)
_NORMALISED_MARKERS: tuple[str, ...] = (NORMALIZED_SUFFIX, "_normalised")


def _normalise_target_export_name(path: Path) -> str:
    """Return a normalised filename for target exports."""

    name = path.name
    lowered_name = name.lower()

    while lowered_name.startswith("."):
        name = name[1:]
        lowered_name = name.lower()

    for suffix in _TEMPORARY_EXPORT_SUFFIXES:
        suffix_lower = suffix.lower()
        while lowered_name.endswith(suffix_lower):
            name = name[: -len(suffix)]
            lowered_name = name.lower()

    for marker in _NORMALISED_MARKERS:
        marker_lower = marker.lower()
        csv_marker = f".csv{marker_lower}"
        marker_csv = f"{marker_lower}.csv"
        if lowered_name.endswith(csv_marker):
            name = name[: -len(marker)]
            lowered_name = name.lower()
            break
        if lowered_name.endswith(marker_csv):
            extension = name[-len(".csv") :]
            name_without_extension = name[: -len(extension)]
            name = f"{name_without_extension[: -len(marker)]}{extension}"
            lowered_name = name.lower()
            break
        if lowered_name.endswith(marker_lower):
            name = name[: -len(marker)]
            lowered_name = name.lower()
            break

    canonical = name.lower()

    if canonical.startswith("output.") and not _matches_expected_target_input_name(
        canonical
    ):
        candidate = canonical[len("output.") :]
        if _matches_expected_target_input_name(candidate):
            canonical = candidate

    return canonical


def _export_stem(name: str) -> str:
    """Return the lowercase stem for ``name`` stripping the CSV extension."""

    lowered = name.lower()
    return lowered[:-4] if lowered.endswith(".csv") else lowered


_UNSUPPORTED_EXPORT_SUFFIXES: tuple[str, ...] = (
    RAW_SUFFIX.lower(),
    "_chembl",
    "_iuphar",
)


def _is_supported_target_export(path: Path) -> bool:
    """Return ``True`` if *path* represents a canonical target export."""

    export_name = _normalise_target_export_name(path)
    stem = _export_stem(export_name)

    if stem.endswith(_UNSUPPORTED_EXPORT_SUFFIXES):
        return False
    if stem == "out" or stem.startswith("out_"):
        return False

    if _matches_expected_target_input_name(export_name):
        return True

    export_lower = export_name.lower()

    if path.suffix.lower() != ".csv":
        if not export_lower.endswith(".csv"):
            return False
        logger.info(
            "target_postprocess_noncanonical_name",
            path=str(path),
            reason="noncanonical_filename",
            canonical=export_name,
        )
        return True

    logger.info(
        "target_postprocess_noncanonical_name",
        path=str(path),
        reason="noncanonical_filename",
        canonical=export_name,
    )
    return True


_OUTPUT_FILENAME_PATTERN = re.compile(r"^output\.(?P<table>.+)_(?P<date>\d{8})$")
_STEM_DATE_PATTERN = re.compile(r"^(?P<table>.+)_(?P<date>\d{8})$")


def _resolve_output_metadata(
    output: Path,
    *,
    date_hint: str | None = None,
    table_hint: str | None = None,
) -> tuple[str, str]:
    """Return ``(table_name, date_tag)`` derived from ``output`` heuristics."""

    candidates = [output.with_suffix("").name, output.stem]
    for candidate in candidates:
        match = _OUTPUT_FILENAME_PATTERN.match(candidate)
        if match:
            return match.group("table"), match.group("date")
        match = _STEM_DATE_PATTERN.match(candidate)
        if match:
            table_value = match.group("table")
            if table_value.startswith("output."):
                table_value = table_value[len("output.") :]
            return table_value, match.group("date")

    sanitized_table = (table_hint or output.stem or "targets").strip() or "targets"
    normalized_table = sanitized_table.replace(" ", "_")
    if date_hint and re.fullmatch(r"\d{8}", date_hint):
        resolved_date = date_hint
    else:
        resolved_date = datetime.now(timezone.utc).strftime("%Y%m%d")
    return normalized_table, resolved_date


def run_target_postprocess_if_requested(
    source: Path,
    *,
    cfg: Config,
    args: argparse.Namespace | None,
    context: IsoformPostprocessContext | None = None,
    ambiguous_classifications: int | None = None,
) -> "PostprocessingPipelineResult" | None:
    """Execute the consolidated target postprocess pipeline when enabled."""

    postprocess_enabled = bool(getattr(args, "postprocess", False)) if args else False
    if not postprocess_enabled:
        logger.info("Postprocessing skipped (flag --postprocess not set)")
        return None

    if not _is_supported_target_export(source):
        logger.info(
            "target_postprocess_skipped",
            path=str(source),
            reason="unsupported_export_name",
        )
        return None

    stem = _export_stem(_normalise_target_export_name(source))
    if stem.endswith("_uniprot"):
        logger.info(
            "target_postprocess_skipped",
            path=str(source),
            reason="unsupported_export_name",
        )
        return None

    destination = source.with_name("output_postprocessed.targets.csv")

    try:
        from library.postprocess import (
            PostprocessingPipelineConfig,
            get_csv_runtime_config,
            get_pipeline_config,
            run_postprocessing_pipeline,
        )
        from library.postprocessing.targets import (
            TARGET_SCHEMA as _TARGET_POSTPROCESS_SCHEMA,
            run_target_pipeline as _run_target_postprocess,
            validate_targets as _validate_target_postprocess,
        )
    except ImportError as exc:  # pragma: no cover - optional dependency missing
        logger.error(
            "target_postprocess_unavailable",
            path=str(source),
            error=str(exc),
        )
        return None

    config_override = getattr(args, "config", None) if args is not None else None
    pipeline_config = get_pipeline_config("targets", config_override)
    params = dict(pipeline_config.params or {})
    runtime_params = dict(params.get("runtime", {}))

    rerun_postprocess = bool(getattr(args, "rerun_postprocess", False)) if args else False
    partial_run = bool(getattr(args, "partial_run", False)) if args else False
    if rerun_postprocess:
        runtime_params["rerun_postprocess"] = True
    if partial_run:
        runtime_params["partial_run"] = True
    if runtime_params:
        params["runtime"] = runtime_params
        pipeline_config = replace(pipeline_config, params=params)

    csv_runtime_cfg = get_csv_runtime_config(pipeline_config)
    runtime_cfg = PostprocessingPipelineConfig(
        pipeline_config=pipeline_config,
        csv_runtime_config=csv_runtime_cfg,
        runner=_run_target_postprocess,
        validator=_validate_target_postprocess,
        schema=_TARGET_POSTPROCESS_SCHEMA,
        logger=logger,
    )

    try:
        result = run_postprocessing_pipeline(
            "targets",
            source,
            destination,
            runtime_cfg,
        )
    except FileNotFoundError as exc:
        logger.error(
            "target_postprocess_failed",
            path=str(source),
            output=str(destination),
            error=str(exc),
        )
        return None
    except Exception as exc:  # pragma: no cover - defensive logging
        logger.exception(
            "target_postprocess_failed",
            path=str(source),
            output=str(destination),
            error=str(exc),
        )
        return None

    metrics = result.metrics
    summary: dict[str, Any] = {}
    pipeline_version_value: str | None = None
    if metrics is not None:
        summary = metrics.summary()
        pipeline_version_value = metrics.pipeline_version

    payload: dict[str, Any] = {
        "path": str(result.output_path),
        "source": str(source),
    }

    if pipeline_version_value:
        payload["pipeline_version"] = pipeline_version_value
    if summary.get("rows") is not None:
        payload["postprocess_rows"] = summary["rows"]
    if summary.get("columns") is not None:
        payload["postprocess_columns"] = summary["columns"]
    if summary.get("duration_s") is not None:
        payload["postprocess_duration_s"] = summary["duration_s"]
    if summary.get("steps") is not None:
        payload["postprocess_steps"] = summary["steps"]

    if result.report_path is not None:
        payload["postprocess_report"] = str(result.report_path)

    logger.info("target_postprocess_done", **payload)
    logger.info("target_organism_postprocess_done", **payload)

    isoform_payload = dict(payload)
    if context and context.http_requests is not None:
        isoform_payload["http_requests"] = context.http_requests
    if ambiguous_classifications is not None:
        isoform_payload["ambiguous_classifications"] = ambiguous_classifications
    logger.info("target_isoform_postprocess_done", **isoform_payload)

    logger.info("target_names_postprocess_done", **payload)
    logger.info("target_iuphar_postprocess_done", **payload)

    return result


def _resolve_parameter(
    namespace: argparse.Namespace,
    cfg_section: Any,
    attr: str,
    *,
    dest: str | None = None,
    default_value: object,
    cli_overrides: Collection[str] | None = None,
    scope: str,
    fallback: tuple[object, str] | None = None,
) -> ParameterLogEntry:
    dest_name = dest or attr
    if cli_overrides is not None and dest_name in cli_overrides:
        value = getattr(namespace, dest_name)
        source = "cli"
    else:
        current = getattr(cfg_section, attr, default_value)
        if current != default_value:
            value = current
            source = "config"
        elif fallback is not None:
            value, source = fallback
        else:
            value = default_value
            source = "default"
    setattr(namespace, dest_name, value)
    if hasattr(cfg_section, attr):
        setattr(cfg_section, attr, value)
    return ParameterLogEntry(scope=scope, name=attr, value=value, source=source)


def _resolve_target_parameters(
    command: str, cfg: Config, args: argparse.Namespace
) -> list[ParameterLogEntry]:
    cli_overrides: Collection[str] = getattr(args, "_cli_overrides", set())
    entries: list[ParameterLogEntry] = []

    def _apply_mode(
        mode: str,
        section: Any,
        *,
        dest_prefix: str = "",
        scope_prefix: str = "",
        fallback_map: dict[str, ParameterLogEntry] | None = None,
    ) -> None:
        defaults: ModeDefaults = TARGET_MODE_DEFAULTS[mode]
        prefix = f"{dest_prefix}_" if dest_prefix else ""
        scope = scope_prefix or mode
        for name in ("column", "chunk_size", "timeout", "limit", "offset"):
            dest = f"{prefix}{name}"
            fallback = None
            if fallback_map is not None and name in fallback_map:
                entry = fallback_map[name]
                fallback = (entry.value, entry.source)
            entry = _resolve_parameter(
                args,
                section,
                name,
                dest=dest,
                default_value=getattr(defaults, name),
                cli_overrides=cli_overrides,
                scope=f"{scope}.{name}" if scope_prefix else scope,
                fallback=fallback,
            )
            entries.append(entry)

    if command in {"chembl", "uniprot", "iuphar"}:
        section = getattr(cfg.target, command)
        _apply_mode(command, section)
    elif command == "all":
        all_section = cfg.target.all
        _apply_mode("all", all_section)
        fallback_lookup = {
            entry.name: entry for entry in entries if entry.scope == "all"
        }
        _apply_mode(
            "chembl",
            cfg.target.chembl,
            dest_prefix="chembl",
            scope_prefix="all.chembl",
            fallback_map=fallback_lookup,
        )
        _apply_mode(
            "uniprot",
            cfg.target.uniprot,
            dest_prefix="uniprot",
            scope_prefix="all.uniprot",
        )
        _apply_mode(
            "iuphar",
            cfg.target.iuphar,
            dest_prefix="iuphar",
            scope_prefix="all.iuphar",
        )

    return entries


_RUSSIAN_KEYBOARD_MAP: dict[str, str] = {
    "q": "й",
    "w": "ц",
    "e": "у",
    "r": "к",
    "t": "е",
    "y": "н",
    "u": "г",
    "i": "ш",
    "o": "щ",
    "p": "з",
    "[": "х",
    "]": "ъ",
    "a": "ф",
    "s": "ы",
    "d": "в",
    "f": "а",
    "g": "п",
    "h": "р",
    "j": "о",
    "k": "л",
    "l": "д",
    ";": "ж",
    "'": "э",
    "z": "я",
    "x": "ч",
    "c": "с",
    "v": "м",
    "b": "и",
    "n": "т",
    "m": "ь",
    ",": "б",
    ".": "ю",
    "`": "ё",
}


def _translate_keyboard_layout(value: str) -> str:
    """Return ``value`` as if typed with the Russian keyboard layout."""

    translated: list[str] = []
    for char in value:
        lower = char.lower()
        mapped = _RUSSIAN_KEYBOARD_MAP.get(lower)
        if mapped is None:
            translated.append(char)
            continue
        translated.append(mapped.upper() if char.isupper() else mapped)
    return "".join(translated)


def _keyboard_aliases(command: str) -> tuple[str, ...]:
    """Return alternative spellings caused by the Russian keyboard layout."""

    base = _translate_keyboard_layout(command)
    variants: list[str] = []
    for candidate in (base, base.capitalize(), base.upper()):
        if candidate != command and candidate not in variants:
            variants.append(candidate)
    return tuple(variants)


COMMAND_KEYBOARD_ALIASES: dict[str, tuple[str, ...]] = {}
for _command_choice in COMMAND_CHOICES:
    _aliases = _keyboard_aliases(_command_choice)
    if _aliases:
        COMMAND_KEYBOARD_ALIASES[_command_choice] = _aliases

COMMAND_ALIAS_TO_CANONICAL: dict[str, str] = {
    alias: command
    for command, aliases in COMMAND_KEYBOARD_ALIASES.items()
    for alias in aliases
}


@dataclass(frozen=True)
class _UniprotCandidate:
    """Container describing a UniProt identifier candidate for a target row."""

    value: str
    source: str
    original_id: str


@dataclass(frozen=True)
class _UniprotQueryPlan:
    """Deterministic mapping of ChEMBL rows to UniProt identifiers."""

    unique_records: list[dict[str, str]]
    row_candidates: list[list[_UniprotCandidate]]
    row_index: list[Any]
    candidate_columns: list[str]


def _split_uniprot_tokens(value: str) -> Iterator[str]:
    """Yield cleaned UniProt identifiers from a pipe-delimited ``value``."""

    for token in value.split("|"):
        token = token.strip()
        if token and any(char.isalnum() for char in token):
            yield token


def _collect_uniprot_candidate_columns(df: pd.DataFrame, cfg: Config) -> list[str]:
    """Return ordered list of columns potentially holding UniProt accessions."""

    primary = cfg.target.all.uniprot_column
    ordered: list[str] = []
    if primary in df.columns:
        ordered.append(primary)

    preferred = ["uniprot_id", "mapping_uniprot_id"]
    for column in preferred:
        if column != primary and column in df.columns and column not in ordered:
            ordered.append(column)

    extra = [
        column
        for column in df.columns
        if column not in ordered
        and any(keyword in column.lower() for keyword in ("uniprot", "accession"))
    ]
    ordered.extend(extra)
    return ordered


def _ensure_merge_column_present(
    df: pd.DataFrame, merge_column: str, cfg: Config
) -> pd.DataFrame:
    """Return ``df`` ensuring that ``merge_column`` exists or fall back to aliases."""

    if merge_column in df.columns:
        return df

    candidate_columns = [
        column
        for column in _collect_uniprot_candidate_columns(df, cfg)
        if column != merge_column and column in df.columns
    ]

    for source in candidate_columns:
        series = df[source]
        cleaned = series.fillna("").astype(str).map(str.strip)
        if cleaned.replace("", pd.NA).dropna().empty:
            continue
        logger.warning(
            "uniprot_merge_column_alias",
            configured=merge_column,
            alias=source,
        )
        aliased = df.copy()
        aliased[merge_column] = cleaned.astype(object)
        return aliased

    logger.error(
        "missing_uniprot_column",
        configured=merge_column,
        aliases=candidate_columns,
        available=sorted(df.columns.tolist()),
    )
    raise PipelineError(
        "Unable to locate a UniProt column. Provide 'uniprot_id' in the input file "
        "or set target.all.uniprot_column to an existing alias such as a column "
        "containing 'uniprot' or 'accession'."
        f" Missing configured column '{merge_column}'."
    )


def _build_uniprot_query_plan(df: pd.DataFrame, cfg: Config) -> _UniprotQueryPlan:
    """Create a deterministic plan for querying UniProt based on ``df``."""

    candidate_columns = _collect_uniprot_candidate_columns(df, cfg)
    unique_records: list[dict[str, str]] = []
    seen: set[str] = set()
    row_candidates: list[list[_UniprotCandidate]] = []
    row_index = list(df.index)

    if not candidate_columns:
        return _UniprotQueryPlan(
            unique_records, [[] for _ in row_index], row_index, candidate_columns
        )

    positions = [df.columns.get_loc(column) for column in candidate_columns]
    for row in df.itertuples(index=False, name=None):
        candidates: list[_UniprotCandidate] = []
        row_seen: set[str] = set()
        for column, pos in zip(candidate_columns, positions, strict=False):
            raw_value = row[pos]
            if not isinstance(raw_value, str) or not raw_value:
                continue
            for token in _split_uniprot_tokens(raw_value):
                if token in row_seen:
                    continue
                row_seen.add(token)
                candidate = _UniprotCandidate(
                    value=token, source=column, original_id=token
                )
                candidates.append(candidate)
                if token not in seen:
                    seen.add(token)
                    unique_records.append(
                        {
                            "uniprot_id": candidate.value,
                            "original_id": candidate.original_id,
                            "source_column": candidate.source,
                        }
                    )
        row_candidates.append(candidates)

    return _UniprotQueryPlan(
        unique_records, row_candidates, row_index, candidate_columns
    )


def _resolve_uniprot_matches(
    plan: _UniprotQueryPlan, uniprot_df: pd.DataFrame
) -> pd.Series:
    """Return preferred UniProt identifier for each ChEMBL row in ``plan``."""

    if not plan.row_index:
        return pd.Series(dtype=object)

    lookup_column = "uniprot_id"
    if lookup_column not in uniprot_df.columns:
        return pd.Series(
            [UNIPROT_MISSING_VALUE for _ in plan.row_index],
            index=plan.row_index,
            dtype=object,
        )

    cleaned = uniprot_df[lookup_column].dropna().astype(str).map(str.strip)

    if "original_id" in uniprot_df.columns:
        original_series = (
            uniprot_df["original_id"].fillna("").astype(str).map(str.strip)
        )
        candidate_map: dict[str, str] = {}
        for canonical, original in zip(cleaned, original_series, strict=False):
            if not canonical:
                continue
            candidate_map.setdefault(canonical, canonical)
            if original:
                for token in _split_uniprot_tokens(original):
                    if token and token not in candidate_map:
                        candidate_map[token] = canonical
    else:
        available = {value for value in cleaned if value}
        candidate_map = {value: value for value in available}
    resolved: list[str] = []
    for candidates in plan.row_candidates:
        match = UNIPROT_MISSING_VALUE
        for candidate in candidates:
            mapped = candidate_map.get(candidate.value)
            if mapped:
                match = mapped
                break
        resolved.append(match)

    return pd.Series(resolved, index=plan.row_index, dtype=object)


def _normalized_output_path(base: Path) -> Path:
    """Return deterministic path for the normalized export derived from ``base``."""

    suffix = base.suffix or ".csv"
    return base.with_name(f"{base.stem}{NORMALIZED_SUFFIX}{suffix}")


def _raw_output_path(base: Path) -> Path:
    """Return default path for the raw payload dump derived from ``base``."""

    suffix = base.suffix or ".csv"
    return base.with_name(f"{base.stem}{RAW_SUFFIX}{suffix}")


def _ensure_parent_directory(path: Path, *, cfg: Config) -> None:
    """Ensure the parent directory for ``path`` exists or raise an error."""

    parent = path.parent
    if parent.exists():
        if not parent.is_dir():
            raise NotADirectoryError(f"{parent} is not a directory")
        return
    if cfg.io.exist_ok:
        parent.mkdir(parents=True, exist_ok=True)
    else:
        raise FileNotFoundError(f"{parent} does not exist")


def _prepare_raw_destination(destination: Path, *, cfg: Config) -> None:
    """Ensure the raw dump destination can be written to safely."""

    _ensure_parent_directory(destination, cfg=cfg)
    if not destination.exists():
        return

    try:
        destination.unlink()
    except PermissionError:
        writable_mode = stat.S_IRUSR | stat.S_IWUSR
        writable_mode |= getattr(stat, "S_IWRITE", 0)
        try:
            os.chmod(destination, writable_mode)
        except OSError as exc:  # pragma: no cover - defensive guard
            raise OSError(f"failed to prepare raw dump destination: {exc}") from exc

        try:
            destination.unlink()
        except OSError as exc:
            raise OSError(f"failed to prepare raw dump destination: {exc}") from exc
    except OSError as exc:  # pragma: no cover - defensive guard
        raise OSError(f"failed to prepare raw dump destination: {exc}") from exc


class _RawDumpStreamWriter:
    """Stream ChEMBL raw payloads to disk without accumulating chunks."""

    _CSV_WRITE_MAX_ATTEMPTS = 5
    _CSV_WRITE_RETRY_BASE_SLEEP_S = 0.1

    def __init__(
        self,
        destination: Path,
        *,
        cfg: Config,
        reindex_columns: bool,
    ) -> None:
        self.destination = destination
        self.cfg = cfg
        self.reindex_columns = reindex_columns
        self._suffix = destination.suffix.lower()
        self._rows_written = 0
        self._columns: list[str] | None = None
        self._frames: list[pd.DataFrame] | None = [] if self._is_parquet else None
        _prepare_raw_destination(destination, cfg=cfg)
        self._destination_opened = False

    @property
    def _is_parquet(self) -> bool:
        return self._suffix in {".parquet", ".pq"}

    def _resolve_columns(self, frame: pd.DataFrame) -> None:
        if self._columns is None:
            columns = list(frame.columns)
            if self.reindex_columns:
                columns = sorted(columns)
            self._columns = columns
            return

        new_columns = [
            column for column in frame.columns if column not in self._columns
        ]
        if not new_columns:
            return

        if self._is_parquet or self._rows_written == 0:
            merged = list(self._columns)
            merged.extend(column for column in new_columns if column not in merged)
            if self.reindex_columns:
                merged = sorted(merged)
            self._columns = merged
            return

        raise OSError("raw_dump_inconsistent_columns")

    def write(self, frame: pd.DataFrame) -> None:
        if frame.empty and self._rows_written == 0 and self._columns is None:
            if self.reindex_columns:
                self._columns = []
            else:
                self._columns = list(frame.columns)
            return

        self._resolve_columns(frame)
        working = frame
        if self._columns is not None and not working.empty:
            working = working.reindex(columns=self._columns)

        if self._is_parquet:
            if self._frames is None:
                self._frames = []
            self._frames.append(working.copy())
        else:
            mode = "w" if self._rows_written == 0 else "a"
            header = self._rows_written == 0
            self._write_csv_with_retry(working, mode=mode, header=header)
        self._rows_written += len(working)
        logger.info(
            "raw_dump_written", rows=self._rows_written, path=str(self.destination)
        )
        return self.destination

    def _write_csv_with_retry(
        self, frame: pd.DataFrame, *, mode: str, header: bool
    ) -> None:
        """Persist ``frame`` to ``destination`` handling transient file locks."""

        last_error: OSError | None = None
        for attempt in range(1, self._CSV_WRITE_MAX_ATTEMPTS + 1):
            try:
                frame.to_csv(
                    str(self.destination),
                    mode=mode,
                    header=header,
                    index=False,
                    sep=self.cfg.io.csv_sep,
                    encoding=self.cfg.io.csv_encoding,
                )
                return
            except OSError as exc:
                should_retry = self._should_retry_on_oserror(exc)
                if not should_retry or attempt == self._CSV_WRITE_MAX_ATTEMPTS:
                    raise
                last_error = exc
                sleep_seconds = self._CSV_WRITE_RETRY_BASE_SLEEP_S * attempt
                logger.warning(
                    "raw_dump_write_retry",
                    attempt=attempt,
                    wait_seconds=sleep_seconds,
                    path=str(self.destination),
                    error=str(exc),
                )
                time.sleep(sleep_seconds)
        if last_error is not None:  # pragma: no cover - defensive
            raise last_error

    @staticmethod
    def _should_retry_on_oserror(exc: OSError) -> bool:
        if isinstance(exc, PermissionError):
            return True
        errno_value = getattr(exc, "errno", None)
        return errno_value in {errno.EACCES, errno.EPERM}

    def finalize(self) -> Path:
        """Flush buffered payloads to ``destination`` and return the path."""

        if self._is_parquet:
            frames = self._frames or []
            if frames:
                combined = pd.concat(frames, ignore_index=True)
            else:
                combined = pd.DataFrame(columns=self._columns or [])
            try:
                combined.to_parquet(self.destination, index=False)
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ValueError(
                    "Parquet export requires optional pyarrow or fastparquet"
                ) from exc
        else:
            if not self.destination.exists():
                empty = pd.DataFrame(columns=self._columns or [])
                empty.to_csv(
                    str(self.destination),
                    index=False,
                    sep=self.cfg.io.csv_sep,
                    encoding=self.cfg.io.csv_encoding,
                )

        return self.destination


def _finalize_raw_dump_writer(
    writer: object,
    *,
    logger: Logger,
    destination: Path,
) -> bool:
    """Finalize ``writer`` if it exposes a ``finalize`` method.

    Parameters
    ----------
    writer:
        Writer instance returned by :class:`_RawDumpStreamWriter`.
    logger:
        Structured logger used by the pipeline.
    destination:
        Raw dump destination path for logging context.

    Returns
    -------
    bool
        ``True`` when the writer either finalizes successfully or does not
        expose a ``finalize`` method. ``False`` indicates that finalization
        failed due to an :class:`OSError`.
    """

    finalize = getattr(writer, "finalize", None)
    if finalize is None or not callable(finalize):
        logger.debug(
            "raw_dump_finalize_missing",
            writer_type=type(writer).__name__,
        )
        return True

    try:
        finalize()
    except OSError as exc:
        logger.error(
            "raw_dump_failed",
            error=str(exc),
            exc_info=exc,
            path=str(destination),
        )
        return False

    return True


def _write_raw_dump(
    df: pd.DataFrame | Iterator[pd.DataFrame],
    destination: Path,
    *,
    cfg: Config,
    reindex_columns: bool,
) -> Path:
    """Persist raw payloads to ``destination`` using streaming writes."""

    writer = _RawDumpStreamWriter(destination, cfg=cfg, reindex_columns=reindex_columns)
    if isinstance(df, pd.DataFrame):
        writer.write(df)
    else:
        for chunk in df:
            writer.write(chunk)
    return writer.finalize()


def _pipe_merge(values: Sequence[str | None]) -> str:
    """Return a ``"|"``-joined string of unique, non-empty tokens.

    Parameters
    ----------
    values:
        Sequence of pipe-delimited strings to merge.

    Returns
    -------
    str
        Sorted, unique tokens separated by ``"|"``. Empty inputs yield
        an empty string.

    """
    tokens: set[str] = set()
    for value in values:
        if isinstance(value, str) and value:
            parts = [p.strip() for p in value.split("|") if p.strip()]
            tokens.update(parts)
    return "|".join(sorted(tokens))


def _first_token(value: str | None) -> str:
    """Return the first token from a pipe-delimited ``value``."""
    if isinstance(value, str) and value:
        return value.split("|")[0]
    return ""


def _prefer_primary(
    primary: pd.Series | None, secondary: pd.Series | None
) -> pd.Series:
    """Return a series prioritising ``primary`` values over ``secondary``.

    The function coalesces two series representing the same logical column
    coming from different data sources. Values from ``primary`` are preferred
    unless they are missing or empty, in which case the corresponding entry
    from ``secondary`` is used. Missing inputs yield an empty object-typed
    series to maintain downstream compatibility.
    """

    if primary is None and secondary is None:
        return pd.Series(dtype=object)

    if primary is None:
        return (
            secondary.astype(object)
            if secondary is not None
            else pd.Series(dtype=object)
        )

    if secondary is None:
        return primary.astype(object)

    primary = primary.astype(object)
    secondary = secondary.astype(object)
    result = primary.copy()
    if len(result) != len(secondary):
        secondary = secondary.reindex(result.index)
    mask = result.isna() | (result == "")
    if mask.any():
        result.loc[mask] = secondary.loc[mask]
    return result


def _save_snapshot(df: pd.DataFrame, base: Path, step: str, cfg: Config) -> Path:
    """Write ``df`` to a uniquely named snapshot CSV file with metadata.

    The file is created alongside ``base`` using the pattern
    ``<base>_<step>_<n>.csv`` where ``n`` increments to avoid overwriting
    existing files. Snapshot exports respect the configured CSV separator and
    encoding, and a ``.meta.yaml`` sidecar is generated containing minimal
    provenance details for reproducibility.

    Parameters
    ----------
    df:
        Data frame to serialise.
    base:
        Base path for the output file. Its stem and suffix determine the
        snapshot file name.
    step:
        Descriptive label inserted into the snapshot file name.
    cfg:
        Application configuration providing CSV export options.

    Returns
    -------
    Path
        Path to the written snapshot file.
    """
    stem = base.stem
    suffix = base.suffix or ".csv"
    index = 1
    while True:
        candidate = base.with_name(f"{stem}_{step}_{index}{suffix}")
        if not candidate.exists():
            work = df.copy()
            if work.columns.empty:
                key_cols: list[str] = []
            else:
                key_cols = [
                    column for column in TARGETS_COLUMN_ORDER if column in work.columns
                ]
                if not key_cols:
                    key_cols = list(work.columns)
            csv_path = write_csv_deterministic(
                work,
                candidate,
                col_order=list(df.columns),
                key_cols=key_cols,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                cfg=None,
            )
            stats: Stats = {
                "rows_total": len(df),
                "rows_kept": len(df),
                "rows_dropped": 0,
                "output_sha256": file_sha256(csv_path),
            }
            write_meta_yaml(
                csv_path=csv_path,
                command=" ".join(sys.argv),
                config=_serialize_paths(cfg.to_dict()),
                inputs={"base": str(base), "step": step},
                stats=stats,
                schema="TargetSnapshot",
            )
            return csv_path
        index += 1


def _build_parser_impl() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create and return the top-level CLI argument parser.

    The command line interface is organised into sub-commands for retrieving
    data from individual sources (UniProt, ChEMBL and IUPHAR) as well as a
    convenience ``all`` command that runs all pipelines and merges their
    outputs.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser populated with every sub-command alongside the logging
        configuration used by :func:`main`.
    """

    root, shared, log_cfg = build_root_parser()

    def _add_output_arguments(
        parser_obj: argparse.ArgumentParser, *, defaults: bool
    ) -> None:
        final_default: object = None if defaults else argparse.SUPPRESS
        raw_default: object = None if defaults else argparse.SUPPRESS
        raw_format_default: object = "csv" if defaults else argparse.SUPPRESS
        id_cols_default: object = None if defaults else argparse.SUPPRESS
        no_reindex_default: object = False if defaults else argparse.SUPPRESS
        normalize_default: object = False if defaults else argparse.SUPPRESS
        option_actions = parser_obj._option_string_actions
        if not any(
            alias in option_actions for alias in ("--final-out", "--out", "--output")
        ):
            parser_obj.add_argument(
                "--final-out",
                "--out",
                "--output",
                dest="final_out",
                type=path_argument,
                default=final_default,
                help=(
                    "Destination for the validated target export "
                    "(default: derived from input name)"
                ),
            )
        parser_obj.add_argument(
            "--raw-out",
            dest="raw_out",
            type=path_argument,
            default=raw_default,
            help=(
                "Location for the raw combined dataset prior to final normalisation"
                " (requires --emit-legacy-artifacts)"
            ),
        )
        parser_obj.add_argument(
            "--raw-format",
            dest="raw_format",
            choices=("csv", "parquet"),
            default=raw_format_default,
            help=(
                "Format used when writing the raw dataset (default: csv). "
                "Requires --emit-legacy-artifacts"
            ),
        )
        parser_obj.add_argument(
            "--id-cols",
            dest="id_cols",
            nargs="+",
            default=id_cols_default,
            help="Identifier columns used for deterministic ordering",
        )
        parser_obj.add_argument(
            "--no-reindex-raw",
            dest="no_reindex_raw",
            action="store_true",
            default=no_reindex_default,
            help="Skip column reindexing when exporting the raw dataset",
        )
        parser_obj.add_argument(
            "--normalize-at-export",
            dest="normalize_at_export",
            action=argparse.BooleanOptionalAction,
            default=normalize_default,
            help=(
                "Apply normalisation immediately before writing the final output. "
                "Use --no-normalize-at-export to keep the raw payload."
            ),
        )
        has_out_alias = any(
            "--out" in action.option_strings for action in parser_obj._actions
        )
        if not has_out_alias:
            parser_obj.add_argument(
                "--out",
                dest="final_out_alias",
                type=path_argument,
                default=argparse.SUPPRESS,
                help=argparse.SUPPRESS,
            )

    _add_output_arguments(root, defaults=True)
    _add_output_arguments(shared, defaults=False)

    root.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser = argparse.ArgumentParser(
        description="Target data utilities",
        parents=[root],
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    def _add_common_cli_options(
        parser_obj: argparse.ArgumentParser,
        *,
        defaults: ModeDefaults,
        column_help: str,
        column_choices: Sequence[str] | None = None,
    ) -> None:
        selection = parser_obj.add_argument_group("Selection")
        selection.add_argument(
            "--column",
            dest="column",
            action=StoreWithSource,
            default=defaults.column,
            choices=column_choices,
            help=f"{column_help} (default: {defaults.column})",
        )
        selection.add_argument(
            "--limit",
            dest="limit",
            action=StoreWithSource,
            default=None,
            type=positive_int,
            help="Maximum number of identifiers to process (omit for no limit)",
        )
        selection.add_argument(
            "--offset",
            dest="offset",
            action=StoreWithSource,
            default=defaults.offset,
            type=_non_negative_int,
            help="Number of identifiers to skip before processing",
        )
        execution = parser_obj.add_argument_group("Execution")
        execution.add_argument(
            "--chunk-size",
            dest="chunk_size",
            action=StoreWithSource,
            default=defaults.chunk_size,
            type=positive_int,
            help=f"Identifiers processed per request (default: {defaults.chunk_size})",
        )
        execution.add_argument(
            "--timeout",
            dest="timeout",
            action=StoreWithSource,
            default=defaults.timeout,
            type=_positive_float,
            help=f"Timeout in seconds for API calls (default: {defaults.timeout})",
        )
        execution.add_argument(
            "--postprocess",
            dest="postprocess",
            action=argparse.BooleanOptionalAction,
            default=False,
            help="Enable target postprocessing after the main pipeline",
        )

    # ----------------------------
    # UniProt sub-command
    # ----------------------------
    uniprot = subparsers.add_parser(
        "uniprot",
        parents=[shared],
        help="Extract information for UniProt accessions",
        aliases=list(COMMAND_KEYBOARD_ALIASES.get("uniprot", ())),
    )
    _add_common_cli_options(
        uniprot,
        defaults=TARGET_MODE_DEFAULTS["uniprot"],
        column_help="Column in the input CSV containing UniProt accessions",
        column_choices=["uniprot_id", "mapping_uniprot_id"],
    )
    uniprot_sources = uniprot.add_argument_group("Data sources")
    uniprot_sources.add_argument(
        "--data-dir",
        dest="data_dir",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config target.uniprot.data_dir)"
        ),
    )
    uniprot.set_defaults(func=run_uniprot)
    uniprot.set_defaults(disable_gtop=False)

    uniprot_network = uniprot.add_argument_group("Network controls")
    uniprot_network.add_argument(
        "--disable-gtop",
        dest="disable_gtop",
        action="store_true",
        help=(
            "Skip Guide-to-Pharmacology enrichment when retrieving UniProt "
            "data to avoid external HTTP requests"
        ),
    )

    # ----------------------------
    # ChEMBL sub-command
    # ----------------------------
    chembl = subparsers.add_parser(
        "chembl",
        parents=[shared],
        help="Retrieve target information from ChEMBL",
        conflict_handler="resolve",
        aliases=list(COMMAND_KEYBOARD_ALIASES.get("chembl", ())),
    )
    chembl.set_defaults(normalize_at_export=True)
    _add_common_cli_options(
        chembl,
        defaults=TARGET_MODE_DEFAULTS["chembl"],
        column_help="Column name in the input CSV containing ChEMBL identifiers",
    )
    chembl.set_defaults(func=run_chembl)

    # ----------------------------
    # IUPHAR sub-command
    # ----------------------------
    iuphar = subparsers.add_parser(
        "iuphar",
        parents=[shared],
        help="Map UniProt accessions to IUPHAR classifications",
        aliases=list(COMMAND_KEYBOARD_ALIASES.get("iuphar", ())),
    )
    _add_common_cli_options(
        iuphar,
        defaults=TARGET_MODE_DEFAULTS["iuphar"],
        column_help="Column in the input CSV containing UniProt accessions",
    )
    iuphar_sources = iuphar.add_argument_group("Data sources")
    iuphar_sources.add_argument(
        "--target-csv",
        dest="target_csv",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config target.iuphar.target_csv)"
        ),
    )
    iuphar_sources.add_argument(
        "--family-csv",
        dest="family_csv",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config target.iuphar.family_csv)"
        ),
    )
    iuphar.set_defaults(func=run_iuphar)

    # ----------------------------
    # Combined pipeline
    # ----------------------------
    all_cmd = subparsers.add_parser(
        "all",
        parents=[shared],
        help="Run ChEMBL, UniProt and IUPHAR pipelines and merge results",
        aliases=list(COMMAND_KEYBOARD_ALIASES.get("all", ())),
    )
    _add_common_cli_options(
        all_cmd,
        defaults=TARGET_MODE_DEFAULTS["all"],
        column_help="Column containing ChEMBL target identifiers",
    )
    outputs_group = all_cmd.add_argument_group("Intermediate outputs")
    outputs_group.add_argument(
        "--chembl-out",
        dest="chembl_out",
        type=path_argument,
        help="Optional path to save intermediate ChEMBL data",
    )
    outputs_group.add_argument(
        "--uniprot-out",
        dest="uniprot_out",
        type=path_argument,
        help="Optional path to save intermediate UniProt data",
    )
    outputs_group.add_argument(
        "--iuphar-out",
        dest="iuphar_out",
        type=path_argument,
        help="Optional path to save intermediate IUPHAR data",
    )
    network_group = all_cmd.add_argument_group("Network controls")
    network_group.add_argument(
        "--disable-gtop",
        dest="disable_gtop",
        action="store_true",
        help=(
            "Skip Guide-to-Pharmacology enrichment during the combined "
            "pipeline to avoid external HTTP requests"
        ),
    )
    all_sources = all_cmd.add_argument_group("Data sources")
    all_sources.add_argument(
        "--data-dir",
        dest="data_dir",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config target.all.data_dir)"
        ),
    )
    all_sources.add_argument(
        "--target-csv",
        dest="target_csv",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config target.all.target_csv)"
        ),
    )
    all_sources.add_argument(
        "--family-csv",
        dest="family_csv",
        type=path_argument,
        default=None,
        action=StoreWithSource,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config target.all.family_csv)"
        ),
    )
    merge_group = all_cmd.add_argument_group("Merge behaviour")
    merge_group.add_argument(
        "--uniprot-column",
        dest="uniprot_column",
        action=StoreWithSource,
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help=(
            "Column from the ChEMBL output used for UniProt processing "
            "(default: config target.all.uniprot_column)"
        ),
    )
    overrides_chembl = all_cmd.add_argument_group("ChEMBL overrides")
    overrides_chembl.add_argument(
        "--chembl-column",
        dest="chembl_column",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["chembl"].column,
        help=(
            "Override the column used for ChEMBL requests "
            f"(default: {TARGET_MODE_DEFAULTS['chembl'].column})"
        ),
    )
    overrides_chembl.add_argument(
        "--chembl-chunk-size",
        dest="chembl_chunk_size",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["chembl"].chunk_size,
        type=positive_int,
        help=(
            "Override the ChEMBL chunk size "
            f"(default: {TARGET_MODE_DEFAULTS['chembl'].chunk_size})"
        ),
    )
    overrides_chembl.add_argument(
        "--chembl-timeout",
        dest="chembl_timeout",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["chembl"].timeout,
        type=_positive_float,
        help=(
            "Override the ChEMBL request timeout "
            f"(default: {TARGET_MODE_DEFAULTS['chembl'].timeout})"
        ),
    )
    overrides_chembl.add_argument(
        "--chembl-limit",
        dest="chembl_limit",
        action=StoreWithSource,
        default=None,
        type=positive_int,
        help="Override the ChEMBL record limit (omit for no limit)",
    )
    overrides_chembl.add_argument(
        "--chembl-offset",
        dest="chembl_offset",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["chembl"].offset,
        type=_non_negative_int,
        help=(
            "Override the number of ChEMBL records to skip "
            f"(default: {TARGET_MODE_DEFAULTS['chembl'].offset})"
        ),
    )
    overrides_uniprot = all_cmd.add_argument_group("UniProt overrides")
    overrides_uniprot.add_argument(
        "--uniprot-chunk-size",
        dest="uniprot_chunk_size",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["uniprot"].chunk_size,
        type=positive_int,
        help=(
            "Override the UniProt chunk size "
            f"(default: {TARGET_MODE_DEFAULTS['uniprot'].chunk_size})"
        ),
    )
    overrides_uniprot.add_argument(
        "--uniprot-timeout",
        dest="uniprot_timeout",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["uniprot"].timeout,
        type=_positive_float,
        help=(
            "Override the UniProt request timeout "
            f"(default: {TARGET_MODE_DEFAULTS['uniprot'].timeout})"
        ),
    )
    overrides_uniprot.add_argument(
        "--uniprot-limit",
        dest="uniprot_limit",
        action=StoreWithSource,
        default=None,
        type=positive_int,
        help="Override the UniProt record limit (omit for no limit)",
    )
    overrides_uniprot.add_argument(
        "--uniprot-offset",
        dest="uniprot_offset",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["uniprot"].offset,
        type=_non_negative_int,
        help=(
            "Override the number of UniProt records to skip "
            f"(default: {TARGET_MODE_DEFAULTS['uniprot'].offset})"
        ),
    )
    overrides_iuphar = all_cmd.add_argument_group("IUPHAR overrides")
    overrides_iuphar.add_argument(
        "--iuphar-column",
        dest="iuphar_column",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["iuphar"].column,
        help=(
            "Override the column used for IUPHAR lookups "
            f"(default: {TARGET_MODE_DEFAULTS['iuphar'].column})"
        ),
    )
    overrides_iuphar.add_argument(
        "--iuphar-chunk-size",
        dest="iuphar_chunk_size",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["iuphar"].chunk_size,
        type=positive_int,
        help=(
            "Override the IUPHAR chunk size "
            f"(default: {TARGET_MODE_DEFAULTS['iuphar'].chunk_size})"
        ),
    )
    overrides_iuphar.add_argument(
        "--iuphar-timeout",
        dest="iuphar_timeout",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["iuphar"].timeout,
        type=_positive_float,
        help=(
            "Override the IUPHAR processing timeout "
            f"(default: {TARGET_MODE_DEFAULTS['iuphar'].timeout})"
        ),
    )
    overrides_iuphar.add_argument(
        "--iuphar-limit",
        dest="iuphar_limit",
        action=StoreWithSource,
        default=None,
        type=positive_int,
        help="Override the IUPHAR record limit (omit for no limit)",
    )
    overrides_iuphar.add_argument(
        "--iuphar-offset",
        dest="iuphar_offset",
        action=StoreWithSource,
        default=TARGET_MODE_DEFAULTS["iuphar"].offset,
        type=_non_negative_int,
        help=(
            "Override the number of IUPHAR records to skip "
            f"(default: {TARGET_MODE_DEFAULTS['iuphar'].offset})"
        ),
    )
    all_cmd.set_defaults(func=run_all)
    all_cmd.set_defaults(disable_gtop=False)

    parser.subparsers_map = {  # type: ignore[attr-defined]
        "uniprot": uniprot,
        "chembl": chembl,
        "iuphar": iuphar,
        "all": all_cmd,
    }

    return parser, log_cfg


def run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``uniprot`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration containing UniProt-specific overrides and
        shared IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser` for the
        ``uniprot`` sub-command.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate failures while reading the
        input CSV, interacting with the UniProt data directory or producing the
        derived artefacts. Input validation errors are logged and converted into
        a failure code.
    """
    disable_gtop_cli = bool(getattr(args, "disable_gtop", False))
    gtop_enabled = (
        cfg.target.uniprot.enable_gtop
        and getattr(cfg.iuphar, "enable", True)
        and not disable_gtop_cli
    )
    if not gtop_enabled:
        if disable_gtop_cli:
            reason = "cli"
        elif not cfg.target.uniprot.enable_gtop:
            reason = "config"
        else:
            reason = "source"
        logger.info("gtop_enrichment_disabled", reason=reason)
        gtop_cfg = cfg.iuphar.model_copy()
        gtop_cfg.enable = False
    else:
        gtop_cfg = cfg.iuphar

    limit = cfg.target.uniprot.limit
    if limit is not None and limit < 1:
        logger.error(
            "invalid_limit",
            section="target.uniprot.limit",
            limit=limit,
        )
        return 1

    output_path: Path | None = None
    export_path: Path | None = None
    resolved_export_path: Path | None = None

    emit_legacy = bool(getattr(args, "emit_legacy_artifacts", False))

    try:
        df = pd.read_csv(
            args.input_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        column = cfg.target.uniprot.column
        if column not in df.columns:
            raise ValueError(f"column '{column}' not found in {args.input_csv}")
        df = df.fillna("")
        df = df[(df[column].str.strip() != "") & (df[column] != "#N/A")].reset_index(
            drop=True
        )
        offset = cfg.target.uniprot.offset
        if offset:
            original_rows = len(df)
            df = df.iloc[offset:].reset_index(drop=True)
            logger.info("process_offset", offset=min(offset, original_rows))
        if limit is not None:
            df = df.head(limit)
            logger.info("process_limit", limit=len(df))
        ids = df[column].to_numpy(copy=False)
        rows_total = len(ids)

        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
        ) as tmp:
            tmp_path = Path(tmp.name)

        pd.DataFrame({"uniprot_id": ids}).to_csv(
            tmp_path,
            index=False,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )

        final_out_attr = getattr(args, "final_out", None)
        if final_out_attr in (None, argparse.SUPPRESS):
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
        else:
            output_path = Path(final_out_attr)
            if not isinstance(final_out_attr, Path):
                args.final_out = output_path
        data_dir = cfg.target.uniprot.data_dir
        uu.init_session(cfg.api, cfg.retry)
        try:
            uu.process(
                input_csv=str(tmp_path),
                output_csv=str(output_path),
                data_dir=data_dir,
                cfg=cfg.uniprot,
                gtop_cfg=gtop_cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        out_df = pd.read_csv(
            output_path, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        rows_kept = len(out_df)
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(
                1,
                "mapping_uniprot_id",
                df["mapping_uniprot_id"].to_numpy(copy=False),
            )
        csv_path = io.write_csv(
            out_df,
            output_path,
            cfg=cfg,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            key_cols=["uniprot_id"],
        )
        export_path = Path(csv_path)
        run_target_postprocess_if_requested(
            export_path,
            cfg=cfg,
            args=args,
            context=IsoformPostprocessContext(args=args),
        )
        rows_dropped = max(rows_total - rows_kept, 0)
        stats: Stats = {
            "rows_total": rows_total,
            "rows_kept": rows_kept,
            "rows_dropped": rows_dropped,
            "output_sha256": file_sha256(csv_path),
        }
        inputs = {"input_csv": str(args.input_csv)}
        if cfg.target.uniprot.data_dir:
            inputs["data_dir"] = str(cfg.target.uniprot.data_dir)
        if emit_legacy:
            write_meta_yaml(
                csv_path=csv_path,
                command=" ".join(sys.argv),
                config=_serialize_paths(cfg.to_dict()),
                inputs=inputs,
                stats=stats,
                schema="UniProtExport",
            )
        else:
            logger.info(
                "legacy_metadata_skipped",
                path=str(csv_path),
                reason="emit_legacy_artifacts_disabled",
            )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "uniprot_processing_failed",
            error=str(exc),
            exc_info=exc,
            input=str(args.input_csv),
            output=str(output_path) if output_path is not None else None,
        )
        return 1
    doc_quality_cfg = cfg.system.doc_quality
    try:
        if is_quality_enabled(doc_quality_cfg):
            if export_path is None:
                raise RuntimeError("export path not available for quality analysis")
            resolved_export_path = export_path.resolve()
            table_quality = build_table_quality_hook(
                doc_quality_cfg,
                table_name=resolved_export_path.with_suffix(""),
                destination=resolved_export_path.parent,
            )
            table_quality(out_df)
    except Exception as exc:
        quality_log_path: Path | None = None
        if resolved_export_path is not None:
            quality_log_path = resolved_export_path
        elif export_path is not None:
            quality_log_path = export_path.resolve()
        elif output_path is not None:
            quality_log_path = output_path.resolve()
        logger.exception(
            "quality_report_failed",
            error=str(exc),
            path=str(quality_log_path) if quality_log_path is not None else None,
            exc=exc,
        )
        return 1
    return 0


def _limited_ids(ids_iter: Iterator[str], limit: int) -> Iterator[str]:
    """Yield up to ``limit`` identifiers while logging the processed count."""

    count = 0
    try:
        for target_id in islice(ids_iter, limit):
            count += 1
            yield target_id
    finally:
        logger.info("process_limit", limit=count)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration providing ChEMBL client settings, retry
        policy and CSV export behaviour.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that reading identifiers,
        fetching data from the ChEMBL API, validating the result or writing the
        CSV artefact failed. All errors are logged with structured context.
    """
    limit = cfg.target.chembl.limit
    if limit is not None and limit < 1:
        logger.error(
            "invalid_limit",
            section="target.chembl.limit",
            limit=limit,
        )
        return 1

    try:
        ids_iter = io.read_ids(
            args.input_csv, column=cfg.target.chembl.column, cfg=cfg.io
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            exc_info=exc,
            path=str(args.input_csv),
        )
        return 1

    offset = cfg.target.chembl.offset
    base_ids_iter: Iterator[str] = ids_iter
    if offset:
        base_ids_iter = islice(base_ids_iter, offset, None)
        logger.info("process_offset", offset=offset)
    if limit is not None:
        base_ids_iter = islice(base_ids_iter, limit)

    processed_ids = 0

    def _counted_ids() -> Iterator[str]:
        nonlocal processed_ids
        for target_id in base_ids_iter:
            processed_ids += 1
            yield target_id

    counted_ids_iter = _counted_ids()

    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        base_output = Path(
            io.default_output_path(
                args.input_csv,
                cfg.io,
                date=getattr(args, "date", None),
            )
        )
        args.final_out = base_output
    else:
        base_output = Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = base_output

    raw_candidate = getattr(args, "raw_out", None)
    if raw_candidate in (None, argparse.SUPPRESS):
        legacy_raw_candidate = getattr(args, "raw_output", None)
        if legacy_raw_candidate not in (None, argparse.SUPPRESS):
            raw_candidate = legacy_raw_candidate
    raw_path_override = (
        Path(raw_candidate) if raw_candidate not in (None, argparse.SUPPRESS) else None
    )
    raw_destination = raw_path_override or _raw_output_path(base_output)
    normalized_output = _normalized_output_path(base_output)

    normalize_flag = getattr(args, "normalize_at_export", None)
    if normalize_flag in (None, argparse.SUPPRESS):
        normalize_at_export = True
    else:
        normalize_at_export = bool(normalize_flag)
    reindex_raw = not bool(getattr(args, "no_reindex_raw", False))

    fetched_rows_total = 0
    raw_dump_rows_total = 0
    chembl_http_requests = 0

    if not normalize_at_export:

        def _raw_fetcher() -> Iterator[pd.DataFrame]:
            nonlocal fetched_rows_total, raw_dump_rows_total, chembl_http_requests
            try:
                with ETLContext(cfg) as context:
                    client = context.chembl_client

                    def _count_attempt() -> None:
                        nonlocal chembl_http_requests
                        chembl_http_requests += 1

                    batch_iter = cl.iter_target_batches_with_retry(
                        counted_ids_iter,
                        cfg=cfg.api,
                        client=client,
                        mapping_cfg=cfg.uniprot_mapping,
                        chunk_size=cfg.target.chembl.chunk_size,
                        timeout=cfg.target.chembl.timeout,
                        retry_cfg=cfg.target.chembl.batch_retry,
                        log=logger,
                        on_attempt=_count_attempt,
                    )
                    for _, raw_chunk, parsed_chunk in batch_iter:
                        raw_dump_rows_total += len(raw_chunk)
                        fetched_rows_total += len(parsed_chunk)
                        if raw_chunk.empty:
                            continue
                        yield raw_chunk
            except (requests.RequestException, ValueError) as exc:
                logger.error(
                    "chembl_fetch_failed",
                    error=str(exc),
                    exc_info=exc,
                    chunk_size=cfg.target.chembl.chunk_size,
                    timeout=cfg.target.chembl.timeout,
                )
                raise PipelineError(str(exc)) from exc

        def _raw_writer(
            chunks: Iterator[pd.DataFrame],
            destination: Path,
            col_order: Sequence[str] | None,
            key_cols: Sequence[str],
        ) -> Path:
            writer = _RawDumpStreamWriter(
                destination, cfg=cfg, reindex_columns=reindex_raw
            )
            for chunk in chunks:
                writer.write(chunk)
            return writer.finalize()

        failure_path = raw_destination.with_name(
            f"{raw_destination.stem}_failure_cases.csv"
        )
        execution = _run_pipeline_with_meta(
            fetcher=_raw_fetcher,
            schema=None,
            schema_name="raw_target_payload",
            validators=[],
            metadata_hooks=[],
            writer=_raw_writer,
            output_path=raw_destination,
            failure_path=failure_path,
            command=" ".join(sys.argv),
            config_snapshot=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(args.input_csv)},
            key_columns=[],
            table_quality=lambda _: None,
            cfg=cfg,
            logger=logger,
        )
        exit_code_attr = getattr(execution, "exit_code", None)
        exit_code = int(exit_code_attr if exit_code_attr is not None else execution)
        dataset_path = getattr(execution, "dataset_path", None) or raw_destination
        raw_destination = Path(dataset_path)
        if exit_code != 0:
            return exit_code

        destination = base_output
        if raw_destination != destination:
            _ensure_parent_directory(destination, cfg=cfg)
            try:
                shutil.copy2(raw_destination, destination)
            except OSError as exc:
                logger.error(
                    "raw_to_final_copy_failed",
                    error=str(exc),
                    exc_info=exc,
                    source=str(raw_destination),
                    destination=str(destination),
                )
                return 1

        run_target_postprocess_if_requested(
            destination,
            cfg=cfg,
            args=args,
            context=IsoformPostprocessContext(
                args=args,
                http_requests=chembl_http_requests,
            ),
        )
        if limit is not None:
            logger.info("process_limit", limit=processed_ids)
        logger.info(
            "chembl_normalization_skipped",
            output=str(destination),
            rows=raw_dump_rows_total,
        )
        return 0

    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        final_output = Path(
            io.default_output_path(
                args.input_csv,
                cfg.io,
                date=getattr(args, "date", None),
            )
        )
        args.final_out = final_output
    else:
        final_output = Path(final_candidate)
        if not isinstance(final_candidate, Path):
            args.final_out = final_output

    if raw_path_override is None:
        raw_output = final_output
    else:
        raw_output = raw_path_override

    raw_format = str(getattr(args, "raw_format", "csv") or "csv").lower()
    if raw_format not in {"csv", "parquet"}:
        logger.warning("unsupported_raw_format", format=raw_format, fallback="csv")
        raw_format = "csv"

    reindex_raw = not bool(getattr(args, "no_reindex_raw", False))
    # ``normalize_at_export`` already coerced above to avoid divergent values.

    id_cols_value = getattr(args, "id_cols", None)
    if id_cols_value in (None, argparse.SUPPRESS):
        key_columns = ["target_chembl_id"]
    elif isinstance(id_cols_value, str):
        key_columns = [id_cols_value]
    else:
        key_columns = list(id_cols_value) or ["target_chembl_id"]
    failure_path = final_output.with_name(f"{final_output.stem}_failure_cases.csv")

    missing_optional_columns: set[str] = set()
    placeholder_replacements = 0
    post_cleanup_rows_total = 0

    def _prepare_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        nonlocal placeholder_replacements, post_cleanup_rows_total
        prepared, _, missing_optional = _prepare_targets_for_schema(frame)
        if missing_optional and not frame.empty:
            placeholder_replacements += len(frame) * len(missing_optional)
        post_cleanup_rows_total += len(prepared)
        missing_optional_columns.update(missing_optional)
        return prepared

    def _normalize_chunk(frame: pd.DataFrame) -> pd.DataFrame:
        nonlocal raw_dump_rows_total
        normalized_chunk = normalize_targets(frame)
        raw_dump_rows_total += len(normalized_chunk)
        return normalized_chunk

    def _validate_chunk(frame: pd.DataFrame) -> ValidationResult:
        try:
            validation = validate_targets(frame, return_result=True)
        except SchemaErrors as exc:
            validated_subset = getattr(exc, "validated_data", frame)
            return ValidationResult(
                validated_subset,
                exc.failure_cases.copy(),
                "TargetsSchema",
            )
        return validation

    raw_dump_writer = _RawDumpStreamWriter(
        raw_destination, cfg=cfg, reindex_columns=reindex_raw
    )

    def fetcher() -> Iterator[pd.DataFrame]:

        nonlocal fetched_rows_total, raw_dump_rows_total, chembl_http_requests

        try:
            with ETLContext(cfg) as context:
                client = context.chembl_client

                def _count_attempt() -> None:
                    nonlocal chembl_http_requests
                    chembl_http_requests += 1

                batch_iter = cl.iter_target_batches_with_retry(
                    counted_ids_iter,
                    cfg=cfg.api,
                    client=client,
                    mapping_cfg=cfg.uniprot_mapping,
                    chunk_size=cfg.target.chembl.chunk_size,
                    timeout=cfg.target.chembl.timeout,
                    retry_cfg=cfg.target.chembl.batch_retry,
                    log=logger,
                    on_attempt=_count_attempt,
                )
                for _, raw_chunk, parsed_chunk in batch_iter:
                    raw_dump_rows_total += len(raw_chunk)
                    try:
                        raw_dump_writer.write(raw_chunk)
                    except OSError as exc:
                        logger.error(
                            "raw_dump_failed",
                            error=str(exc),
                            exc_info=exc,
                            path=str(raw_destination),
                        )
                        raise PipelineError(str(exc)) from exc
                    if parsed_chunk.empty:
                        continue
                    fetched_rows_total += len(parsed_chunk)
                    yield parsed_chunk
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_fetch_failed",
                error=str(exc),
                exc_info=exc,
                chunk_size=cfg.target.chembl.chunk_size,
                timeout=cfg.target.chembl.timeout,
            )
            raise PipelineError(str(exc)) from exc

    def writer(
        chunks: Iterator[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str] | None,
        key_cols: Sequence[str],
    ) -> Path:
        resolved_keys = list(key_cols) if key_cols else key_columns
        column_order: Sequence[str] | None = col_order if reindex_raw else None

        if raw_format == "csv" and not normalize_at_export:
            raw_path = io.write_csv(
                chunks,
                raw_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=column_order,
                chunksize=cfg.io.csv_chunksize,
            )
            if final_output != raw_output:
                shutil.copy2(raw_path, final_output)
                final_path = final_output
            else:
                final_path = raw_path
            run_target_postprocess_if_requested(
                final_path,
                cfg=cfg,
                args=args,
                context=IsoformPostprocessContext(
                    args=args,
                    http_requests=chembl_http_requests,
                ),
            )
            return final_path

        frames: list[pd.DataFrame] = []
        for chunk in chunks:
            working = chunk
            if column_order is not None:
                working = working.reindex(columns=column_order)
            frames.append(working)

        if frames:
            combined = pd.concat(frames, ignore_index=True)
        else:
            combined = pd.DataFrame(columns=column_order or [])

        if resolved_keys:
            missing_keys = [col for col in resolved_keys if col not in combined.columns]
            if not missing_keys and resolved_keys:
                combined = combined.sort_values(by=list(resolved_keys)).reset_index(
                    drop=True
                )

        if raw_format == "parquet":
            try:
                combined.to_parquet(raw_output, index=False)
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ValueError(
                    "Parquet export requires optional pyarrow or fastparquet"
                ) from exc
            raw_path = raw_output
        else:
            raw_path = io.write_csv(
                combined,
                raw_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=column_order,
            )

        final_df = combined
        if normalize_at_export:
            final_df = normalize_targets(final_df)
            final_df, _, missing_optional = _prepare_targets_for_schema(final_df)
            if missing_optional is not None:
                missing_optional_columns.update(missing_optional)
                missing_optional_columns.update(missing_optional)

        if (
            final_output == raw_output
            and not normalize_at_export
            and raw_format == "parquet"
        ):
            final_path = raw_path
        else:
            final_path = io.write_csv(
                final_df,
                final_output,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=resolved_keys or None,
                col_order=TARGETS_COLUMN_ORDER,
            )
        run_target_postprocess_if_requested(
            final_path,
            cfg=cfg,
            args=args,
            context=IsoformPostprocessContext(
                args=args,
                http_requests=chembl_http_requests,
            ),
        )
        return final_path

    failure_path = normalized_output.with_name(
        f"{normalized_output.stem}_failure_cases.csv"
    )
    doc_quality_cfg = cfg.system.doc_quality
    table_quality = build_table_quality_hook(
        doc_quality_cfg,
        table_name=final_output.with_suffix(""),
        destination=final_output.parent,
    )

    metadata_hooks = [add_pipeline_metadata, _prepare_chunk]
    if not normalize_at_export:
        metadata_hooks.insert(0, normalize_targets)

    execution = _run_pipeline_with_meta(
        fetcher=fetcher,
        schema=TargetsSchema,
        schema_name="TargetsSchema",
        validators=[_validate_chunk],
        metadata_hooks=metadata_hooks,
        writer=writer,
        output_path=raw_output,
        failure_path=failure_path,
        command=" ".join(sys.argv),
        config_snapshot=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        key_columns=key_columns,
        table_quality=table_quality,
        cfg=cfg,
        logger=logger,
        dictionary_resources=(
            "target_uniprot_cache",
            "target_iuphar_target",
            "target_iuphar_family",
        ),
    )
    exit_code_attr = getattr(execution, "exit_code", None)
    exit_code = int(exit_code_attr if exit_code_attr is not None else execution)
    dataset_path = getattr(execution, "dataset_path", None) or raw_output
    raw_output = Path(dataset_path)

    if not _finalize_raw_dump_writer(
        raw_dump_writer,
        logger=logger,
        destination=raw_destination,
    ):
        return 1

    if limit is not None:
        logger.info("process_limit", limit=processed_ids)

    if missing_optional_columns:
        logger.debug(
            "schema_optional_columns_missing",
            columns=sorted(missing_optional_columns),
        )

    logger.info("chembl_stage_rows", stage="fetch", rows=fetched_rows_total)
    logger.info("chembl_stage_rows", stage="raw_dump", rows=raw_dump_rows_total)
    logger.info(
        "chembl_stage_rows",
        stage="post_cleanup",
        rows=post_cleanup_rows_total,
    )
    logger.info(
        "chembl_placeholder_replacements",
        total=placeholder_replacements,
    )

    return exit_code


def run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``iuphar`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration containing paths to the IUPHAR export files
        and shared IO settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate that the IUPHAR sources
        could not be read, the combined dataset failed validation or the CSV
        output could not be written.
    """
    limit = cfg.target.iuphar.limit
    if limit is not None and limit < 1:
        logger.error(
            "invalid_limit",
            section="target.iuphar.limit",
            limit=limit,
        )
        return 1

    tmp_path: Path | None = None
    source_csv = args.input_csv

    output_path_candidate = getattr(args, "output_csv", None)
    if output_path_candidate not in (None, argparse.SUPPRESS):
        output_path = Path(output_path_candidate)
        if not isinstance(output_path_candidate, Path):
            args.output_csv = output_path
        args.final_out = output_path
    else:
        final_out_attr = getattr(args, "final_out", None)
        if final_out_attr in (None, argparse.SUPPRESS):
            output_path = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
            args.final_out = output_path
        else:
            output_path = Path(final_out_attr)
            if not isinstance(final_out_attr, Path):
                args.final_out = output_path

    try:
        df_to_process: pd.DataFrame | None = None
        offset = cfg.target.iuphar.offset
        needs_mutation = (
            limit is not None or offset or cfg.target.iuphar.column != "uniprot_id"
        )
        if needs_mutation:
            df_full = pd.read_csv(
                args.input_csv,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                dtype=str,
            ).fillna("")
            source_column = cfg.target.iuphar.column
            if source_column not in df_full.columns:
                raise ValueError(
                    f"column '{source_column}' not found in {args.input_csv}"
                )
            if source_column != "uniprot_id":
                df_full["uniprot_id"] = df_full[source_column]
            if offset:
                total_rows = len(df_full)
                df_full = df_full.iloc[offset:].reset_index(drop=True)
                logger.info("process_offset", offset=min(offset, total_rows))
            if limit is not None:
                df_full = df_full.head(limit)
                logger.info("process_limit", limit=len(df_full))
            df_to_process = df_full
        if df_to_process is not None:
            from tempfile import NamedTemporaryFile

            with NamedTemporaryFile(delete=False) as tmp:
                tmp_path = Path(tmp.name)
            df_to_process.to_csv(
                tmp_path,
                index=False,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
            source_csv = tmp_path

        data = ii.IUPHARData.from_files(
            target_path=cfg.target.iuphar.target_csv,
            family_path=cfg.target.iuphar.family_csv,
            encoding=cfg.io.csv_encoding,
        )
        data.map_uniprot_file(
            input_path=source_csv,
            output_path=str(output_path),
            encoding=cfg.io.csv_encoding,
            sep=cfg.io.csv_sep,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error(
            "iuphar_processing_failed",
            error=str(exc),
            exc_info=exc,
            input=str(source_csv),
            target_csv=str(cfg.target.iuphar.target_csv),
            family_csv=str(cfg.target.iuphar.family_csv),
        )
        return 1
    finally:
        if tmp_path is not None:
            tmp_path.unlink(missing_ok=True)
    doc_quality_cfg = cfg.system.doc_quality
    table_quality = build_table_quality_hook(
        doc_quality_cfg,
        table_name=output_path.with_suffix(""),
        destination=output_path.parent,
    )
    try:
        table_quality(output_path)
    except Exception as exc:
        logger.exception(
            "quality_report_failed",
            error=str(exc),
            path=str(output_path),
            exc=exc,
        )
        return 1
    return 0


def fetch_chembl(
    cfg: Config,
    input_csv: Path,
    final_out: Path,
    limit: int | None = None,
    *,
    raw_out: Path | None = None,
    raw_format: str = "csv",
    id_cols: Sequence[str] | None = None,
    chunk_size: int | None = None,
    offset: int = 0,
    normalize_at_export: bool = True,
    no_reindex_raw: bool = False,
) -> pd.DataFrame:
    """Fetch target information from ChEMBL.

    Parameters
    ----------
    cfg : Config
        Application configuration used to drive the ChEMBL pipeline.
    input_csv : pathlib.Path
        Source CSV containing target identifiers.
    final_out : pathlib.Path
        Destination path used by :func:`run_chembl` to persist results.
    limit : int, optional
        Maximum number of identifiers to process. ``None`` processes all rows.
    chunk_size : int, optional
        Temporary override for the batch size used when calling the API.
    offset : int, optional
        Number of identifiers to skip before starting the retrieval.

    Returns
    -------
    pandas.DataFrame
        Retrieved ChEMBL data loaded from ``final_out``. The export is
        normalised by default to ensure that downstream steps receive the
        canonical schema (including columns such as ``uniprot_id``).

    Raises
    ------
    RuntimeError
        Raised when :func:`run_chembl` reports a non-zero exit code.
    """

    logger.info("fetch_chembl_start", input=str(input_csv), output=str(final_out))
    chembl_args = argparse.Namespace(
        input_csv=input_csv,
        final_out=final_out,
        output_csv=final_out,
        raw_out=raw_out,
        raw_format=raw_format,
        id_cols=list(id_cols) if id_cols is not None else None,
        offset=offset,
        normalize_at_export=normalize_at_export,
        no_reindex_raw=no_reindex_raw,
    )
    original_limit = cfg.target.chembl.limit
    original_chunk_size = cfg.target.chembl.chunk_size
    if limit is not None:
        cfg.target.chembl.limit = limit
    if chunk_size is not None:
        cfg.target.chembl.chunk_size = chunk_size
    try:
        if run_chembl(cfg, chembl_args) != 0:
            raise RuntimeError("ChEMBL retrieval failed")
    finally:
        if limit is not None:
            cfg.target.chembl.limit = original_limit
        if chunk_size is not None:
            cfg.target.chembl.chunk_size = original_chunk_size
    _normalized_output_path(final_out)
    df = pd.read_csv(
        final_out, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    logger.info("fetch_chembl_done", rows=len(df))
    return df


def fetch_uniprot(
    cfg: Config, chembl_df: pd.DataFrame, output_csv: Path
) -> pd.DataFrame:
    """Retrieve UniProt annotations for targets in ``chembl_df``.

    Parameters
    ----------
    cfg : Config
        Application configuration used to invoke :func:`run_uniprot`.
    chembl_df : pandas.DataFrame
        DataFrame containing ChEMBL target records with at least one UniProt
        identifier column.
    output_csv : pathlib.Path
        Destination path populated by the UniProt pipeline.

    Returns
    -------
    pandas.DataFrame
        UniProt records with an additional ``original_id`` column preserving
        the queried accessions.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_uniprot` returns a non-zero exit code.
    """

    logger.info("fetch_uniprot_start", output=str(output_csv))
    plan = _build_uniprot_query_plan(chembl_df, cfg)
    logger.info(
        "fetch_uniprot_plan",
        rows=len(plan.row_index),
        unique=len(plan.unique_records),
        candidate_columns=plan.candidate_columns,
    )
    if not plan.unique_records:
        logger.info(
            "fetch_uniprot_no_candidates",
            output=str(output_csv),
            rows=len(plan.row_index),
            candidate_columns=plan.candidate_columns,
        )
        empty_df = pd.DataFrame(
            {
                "uniprot_id": pd.Series(dtype=object),
                "original_id": pd.Series(dtype=object),
                "source_column": pd.Series(dtype=object),
                "mapping_uniprot_id": pd.Series(dtype=object),
            }
        )
        write_csv_deterministic(
            empty_df,
            output_csv,
            col_order=[
                "uniprot_id",
                "original_id",
                "source_column",
                "mapping_uniprot_id",
            ],
            key_cols=["uniprot_id"],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        return empty_df

    id_df = pd.DataFrame(plan.unique_records, dtype=object)

    id_df = id_df.copy()
    id_df["__query_order"] = range(len(id_df))
    query_input_df = id_df.drop(columns=["__query_order"], errors="ignore")

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(
        "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
    ) as tmp:
        tmp_path = Path(tmp.name)

    query_input_df.to_csv(
        tmp_path,
        index=False,
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
    )

    uniprot_args = argparse.Namespace(
        input_csv=tmp_path,
        final_out=output_csv,
        output_csv=output_csv,
    )
    orig_dir = cfg.target.uniprot.data_dir
    cfg.target.uniprot.data_dir = cfg.target.all.data_dir
    try:
        if run_uniprot(cfg, uniprot_args) != 0:
            raise RuntimeError("UniProt retrieval failed")
    finally:
        cfg.target.uniprot.data_dir = orig_dir
        tmp_path.unlink(missing_ok=True)

    if not output_csv.exists():
        logger.warning(
            "fetch_uniprot_output_missing",
            path=str(output_csv),
        )
        placeholder = pd.DataFrame(
            {
                "uniprot_id": pd.Series(dtype=object),
                "original_id": pd.Series(dtype=object),
                "source_column": pd.Series(dtype=object),
                "mapping_uniprot_id": pd.Series(dtype=object),
            }
        )
        write_csv_deterministic(
            placeholder,
            output_csv,
            col_order=[
                "uniprot_id",
                "original_id",
                "source_column",
                "mapping_uniprot_id",
            ],
            key_cols=["uniprot_id"],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )

    fetched_df = pd.read_csv(
        output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
    )
    requested_ids = {
        str(record.get("uniprot_id", "")).strip()
        for record in plan.unique_records
        if str(record.get("uniprot_id", "")).strip()
    }
    if "uniprot_id" in fetched_df.columns:
        retrieved_ids = {
            str(value).strip()
            for value in fetched_df["uniprot_id"].dropna()
            if str(value).strip()
        }
    else:
        retrieved_ids = set()
    missing_ids = requested_ids - retrieved_ids
    logger.info(
        "fetch_uniprot_coverage",
        requested=len(requested_ids),
        retrieved=len(retrieved_ids),
        missing=len(missing_ids),
    )
    if missing_ids:
        sample = sorted(missing_ids)[:10]
        logger.debug(
            "fetch_uniprot_missing_ids",
            sample=sample,
            sample_size=len(sample),
            total_missing=len(missing_ids),
        )
    if "uniprot_id" in fetched_df.columns:
        df = id_df.merge(fetched_df, on="uniprot_id", how="left", sort=False)
    else:
        df = id_df.copy()

    if "__query_order" in df.columns:
        df = df.sort_values("__query_order").drop(
            columns=["__query_order"], errors="ignore"
        )
        df = df.reset_index(drop=True)

    left_original = df.pop("original_id_x") if "original_id_x" in df.columns else None
    right_original = df.pop("original_id_y") if "original_id_y" in df.columns else None

    if left_original is not None or right_original is not None:
        df["original_id"] = _prefer_primary(left_original, right_original)
    elif "original_id" not in df.columns:
        df["original_id"] = pd.Series(dtype=object)

    source_column = pd.Series(dtype=object)
    if "source_column_x" in df.columns or "source_column_y" in df.columns:
        left_source = (
            df.pop("source_column_x") if "source_column_x" in df.columns else None
        )
        right_source = (
            df.pop("source_column_y") if "source_column_y" in df.columns else None
        )
        source_column = _prefer_primary(left_source, right_source)
    elif "source_column" in df.columns:
        source_column = df["source_column"].astype(object)
    else:
        source_column = pd.Series(
            [UNIPROT_MISSING_VALUE] * len(df), index=df.index, dtype=object
        )

    if "source_column" not in df.columns:
        df["source_column"] = source_column.reindex(
            df.index, fill_value=UNIPROT_MISSING_VALUE
        )

    if "mapping_uniprot_id" not in df.columns:
        df["mapping_uniprot_id"] = pd.Series(
            UNIPROT_MISSING_VALUE, index=df.index, dtype=object
        )

    mapping_mask = df["source_column"].astype(str).eq("mapping_uniprot_id")
    if mapping_mask.any():
        df.loc[mapping_mask, "mapping_uniprot_id"] = (
            df.loc[mapping_mask, "original_id"]
            .fillna(UNIPROT_MISSING_VALUE)
            .astype(str)
        )
    df["mapping_uniprot_id"] = df["mapping_uniprot_id"].fillna("").astype(str)
    logger.info("fetch_uniprot_done", rows=len(df))
    return df


def _prepare_iuphar_merge_frames(
    combined_df: pd.DataFrame, iuphar_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return frames ready for merging IUPHAR classifications."""

    sanitised_combined = combined_df.drop(
        columns=[col for col in _IUPHAR_OVERRIDE_COLUMNS if col in combined_df.columns],
        errors="ignore",
    )
    existing_cols = set(sanitised_combined.columns)
    keep_columns: list[str] = ["uniprot_id"]
    for column in iuphar_df.columns:
        if column == "uniprot_id":
            continue
        if column in _IUPHAR_OVERRIDE_COLUMNS or column not in existing_cols:
            keep_columns.append(column)
    trimmed_iuphar = iuphar_df.loc[:, keep_columns].copy()
    return sanitised_combined, trimmed_iuphar


def fetch_iuphar(
    cfg: Config,
    chembl_df: pd.DataFrame,
    uniprot_df: pd.DataFrame,
    output_csv: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retrieve IUPHAR classifications for combined target data.

    Parameters
    ----------
    cfg : Config
        Application configuration containing IUPHAR resource paths.
    chembl_df : pandas.DataFrame
        ChEMBL target data obtained from :func:`fetch_chembl`.
    uniprot_df : pandas.DataFrame
        UniProt annotations with an ``original_id`` column provided by
        :func:`fetch_uniprot`.
    output_csv : pathlib.Path
        Destination where :func:`run_iuphar` persists the retrieved
        classifications.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Two data frames: the merged ChEMBL/UniProt input and the IUPHAR
        classification results.

    Raises
    ------
    RuntimeError
        Raised when :func:`run_iuphar` reports a non-zero exit code.
    """

    logger.info("fetch_iuphar_start", output=str(output_csv))
    merge_column = cfg.target.all.uniprot_column
    chembl_df = _ensure_merge_column_present(chembl_df, merge_column, cfg)
    if merge_column not in chembl_df.columns:
        logger.error(
            "invalid_uniprot_column",
            configured=merge_column,
            available=sorted(chembl_df.columns.tolist()),
        )
        raise PipelineError(
            f"Configured UniProt column '{merge_column}' is not present in the ChEMBL target data."
        )

    plan = _build_uniprot_query_plan(chembl_df, cfg)
    lookup_series = _resolve_uniprot_matches(plan, uniprot_df)
    if lookup_series.empty and len(chembl_df.index):
        lookup_series = pd.Series(
            [UNIPROT_MISSING_VALUE] * len(chembl_df.index),
            index=chembl_df.index,
            dtype=object,
        )
    else:
        lookup_series = lookup_series.reindex(
            chembl_df.index, fill_value=UNIPROT_MISSING_VALUE
        )
    lookup_column = "__uniprot_lookup_id"

    if merge_column != "uniprot_id":
        chembl_for_merge = chembl_df.drop(columns=["uniprot_id"], errors="ignore")
    else:
        chembl_for_merge = chembl_df.copy()

    mapping_series = chembl_for_merge.get("mapping_uniprot_id")
    if mapping_series is not None:
        mapping_series = mapping_series.fillna("").astype(str).map(str.strip)
    else:
        mapping_series = pd.Series(
            UNIPROT_MISSING_VALUE, index=chembl_for_merge.index, dtype=object
        )

    resolved_series = lookup_series.reindex(
        chembl_for_merge.index, fill_value=UNIPROT_MISSING_VALUE
    )
    resolved_series = _prefer_primary(resolved_series, mapping_series)
    chembl_for_merge[lookup_column] = resolved_series

    mapping_join_column = "__uniprot_mapping_lookup_id"
    chembl_for_merge[mapping_join_column] = mapping_series.reindex(
        chembl_for_merge.index, fill_value=UNIPROT_MISSING_VALUE
    )

    uniprot_df = uniprot_df.copy()
    if "mapping_uniprot_id" in uniprot_df.columns:
        uniprot_df["mapping_uniprot_id"] = (
            uniprot_df["mapping_uniprot_id"].fillna("").astype(str).map(str.strip)
        )
    else:
        uniprot_df["mapping_uniprot_id"] = pd.Series(
            UNIPROT_MISSING_VALUE, index=uniprot_df.index, dtype=object
        )

    combined_df = pd.merge(
        chembl_for_merge,
        uniprot_df,
        left_on=[lookup_column, mapping_join_column],
        right_on=["uniprot_id", "mapping_uniprot_id"],
        how="left",
        suffixes=("_chembl", "_uniprot"),
    )
    combined_df = combined_df.copy()

    merge_suffix_priority = ("_uniprot", "_chembl")

    def _base_name(column: str) -> str | None:
        for suffix in merge_suffix_priority:
            if column.endswith(suffix):
                return column[: -len(suffix)]
        return None

    coalesced_updates: dict[str, pd.Series] = {}

    def _coalesce_column(df: pd.DataFrame, column: str) -> None:
        preferred: pd.Series | None = None
        for suffix in merge_suffix_priority:
            candidate_name = f"{column}{suffix}"
            if candidate_name not in df.columns:
                continue
            candidate = df.pop(candidate_name).astype(object)
            if preferred is None:
                preferred = candidate
            else:
                preferred = _prefer_primary(preferred, candidate)
        if preferred is not None:
            coalesced_updates[column] = preferred.reindex(df.index)

    critical_columns = set(TARGETS_COLUMN_ORDER)
    critical_columns.update(
        {
            "gene",
            "gene_name",
            "synonyms",
            "component_description",
            "chembl_alternative_name",
            "names",
            "mapping_uniprot_id",
            "ec_numbers",
            "reaction_ec_numbers",
            "ec_code",
        }
    )

    suffixed_columns = [
        column for column in combined_df.columns if _base_name(column) is not None
    ]
    base_columns = {
        _base_name(column) for column in suffixed_columns if _base_name(column)
    }

    for column in sorted(base_columns & critical_columns):
        _coalesce_column(combined_df, column)

    for column in sorted(base_columns - critical_columns):
        _coalesce_column(combined_df, column)

    if coalesced_updates:
        combined_df = combined_df.assign(**coalesced_updates)
        coalesced_updates.clear()

    remaining_suffixes = [
        column for column in combined_df.columns if _base_name(column) is not None
    ]
    if remaining_suffixes:
        rename_map = {
            column: _base_name(column)
            for column in remaining_suffixes
            if _base_name(column)
        }
        combined_df = combined_df.rename(columns=rename_map)

    combined_df = combined_df.drop(
        columns=[
            col for col in [lookup_column, mapping_join_column] if col is not None
        ],
        errors="ignore",
    )
    if "original_id" in combined_df.columns:
        combined_df = combined_df.drop(columns=["original_id"])

    overlap_columns = sorted(
        col
        for col in (
            set(chembl_for_merge.columns) & (set(uniprot_df.columns) - {"original_id"})
        )
        if col is not None
    )
    overlap_updates: dict[str, pd.Series] = {}

    for column in overlap_columns:
        chembl_col = f"{column}_chembl"
        uniprot_col = f"{column}_uniprot"
        chembl_series = (
            combined_df.pop(chembl_col) if chembl_col in combined_df.columns else None
        )
        uniprot_series = (
            combined_df.pop(uniprot_col) if uniprot_col in combined_df.columns else None
        )
        if chembl_series is None and uniprot_series is None:
            continue
        if column == "reaction_ec_numbers":
            chembl_values = (
                chembl_series.reindex(combined_df.index)
                if chembl_series is not None
                else pd.Series(index=combined_df.index, dtype=object)
            )
            uniprot_values = (
                uniprot_series.reindex(combined_df.index)
                if uniprot_series is not None
                else pd.Series(index=combined_df.index, dtype=object)
            )
            overlap_updates[column] = pd.Series(
                [
                    normalize_reaction_ec_numbers([u, c])
                    for u, c in zip(uniprot_values, chembl_values, strict=False)
                ],
                index=combined_df.index,
                dtype=object,
            )
        else:
            overlap_updates[column] = _prefer_primary(
                uniprot_series, chembl_series
            ).reindex(combined_df.index)

    if overlap_updates:
        combined_df = combined_df.assign(**overlap_updates)

    column_updates: dict[str, pd.Series] = {}
    index = combined_df.index

    if "gene" in combined_df.columns:
        gene_series = combined_df["gene"].astype(object)
    else:
        gene_series = pd.Series(
            UNIPROT_MISSING_VALUE,
            index=index,
            dtype=object,
        )
        column_updates["gene"] = gene_series

    ec_number_columns = [
        column for column in combined_df.columns if column.startswith("ec_numbers")
    ]
    ec_numbers_series: pd.Series | None = None
    if ec_number_columns:
        ec_numbers_series = combined_df.apply(
            lambda r: _pipe_merge([r.get(column) for column in ec_number_columns]),
            axis=1,
        )
        column_updates["ec_numbers"] = ec_numbers_series
        extra_ec_columns = [
            column for column in ec_number_columns if column != "ec_numbers"
        ]
        if extra_ec_columns:
            combined_df = combined_df.drop(columns=extra_ec_columns, errors="ignore")

    reaction_ec_columns = [
        column
        for column in combined_df.columns
        if column.startswith("reaction_ec_numbers")
    ]
    reaction_ec_series: pd.Series | None = None
    if reaction_ec_columns:
        reaction_ec_series = combined_df.apply(
            lambda r: normalize_reaction_ec_numbers(
                [r.get(column) for column in reaction_ec_columns]
            ),
            axis=1,
        )
        column_updates["reaction_ec_numbers"] = reaction_ec_series
        extra_reaction_columns = [
            column for column in reaction_ec_columns if column != "reaction_ec_numbers"
        ]
        if extra_reaction_columns:
            combined_df = combined_df.drop(
                columns=extra_reaction_columns, errors="ignore"
            )
    elif "reaction_ec_numbers" in combined_df.columns:
        reaction_ec_series = combined_df["reaction_ec_numbers"].astype(object)

    def _series_or_empty(name: str) -> pd.Series:
        if name in combined_df.columns:
            return combined_df[name].reindex(index).astype(object)
        return pd.Series(index=index, dtype=object)

    synonyms_components = pd.concat(
        [
            _series_or_empty("pref_name"),
            _series_or_empty("component_description"),
            gene_series.reindex(index),
            _series_or_empty("chembl_alternative_name"),
            _series_or_empty("names"),
            _series_or_empty("secondaryAccessionNames"),
        ],
        axis=1,
    )
    column_updates["synonyms"] = synonyms_components.apply(
        lambda r: _pipe_merge(r.tolist()),
        axis=1,
    )

    if ec_numbers_series is None and "ec_numbers" in combined_df.columns:
        ec_numbers_series = combined_df["ec_numbers"].reindex(index).astype(object)
    if ec_numbers_series is None:
        ec_numbers_series = pd.Series(index=index, dtype=object)

    if reaction_ec_series is None:
        reaction_ec_series = pd.Series(index=index, dtype=object)
    else:
        reaction_ec_series = reaction_ec_series.reindex(index).astype(object)

    ec_numbers_values = ec_numbers_series.reindex(index, fill_value=None)
    reaction_values = reaction_ec_series.reindex(index, fill_value=None)
    column_updates["ec_number"] = pd.Series(
        (
            _pipe_merge([ec_val, reaction_val])
            for ec_val, reaction_val in zip(
                ec_numbers_values, reaction_values, strict=False
            )
        ),
        index=index,
        dtype=object,
    )

    column_updates["gene_name"] = gene_series.reindex(index).apply(_first_token)

    if "mapping_uniprot_id" in combined_df.columns:
        column_updates["mapping_uniprot_id"] = (
            combined_df["mapping_uniprot_id"].fillna("").astype(str)
        )

    if column_updates:
        combined_df = combined_df.assign(**column_updates)

    combined_df = combined_df.drop(columns=["ec_numbers"], errors="ignore")

    resolved_for_combined = resolved_series.reindex(
        combined_df.index, fill_value=UNIPROT_MISSING_VALUE
    )
    resolved_for_combined = resolved_for_combined.astype(object)

    if merge_column in combined_df.columns:
        existing_merge = combined_df[merge_column].astype(object)
    else:
        logger.error(
            "missing_configured_uniprot_column",
            configured=merge_column,
            available=sorted(combined_df.columns.tolist()),
        )
        raise PipelineError(
            "Configured UniProt column "
            f"'{merge_column}' is not present in the merged target data."
        )

    if merge_column == "uniprot_id":
        if "mapping_uniprot_id" in combined_df.columns:
            mapping_merge = combined_df["mapping_uniprot_id"].astype(object)
        else:
            mapping_merge = pd.Series(
                UNIPROT_MISSING_VALUE, index=combined_df.index, dtype=object
            )
        preferred_ids = _prefer_primary(resolved_for_combined, mapping_merge)
        preferred_ids = _prefer_primary(preferred_ids, existing_merge)
        combined_df[merge_column] = preferred_ids.fillna(UNIPROT_MISSING_VALUE)
    else:
        preferred_ids = _prefer_primary(existing_merge, resolved_for_combined)
        combined_df[merge_column] = preferred_ids.fillna(UNIPROT_MISSING_VALUE)

    chembl_reactions = chembl_df.get(
        "reaction_ec_numbers", pd.Series(dtype=object)
    ).reindex(combined_df.index, fill_value=UNIPROT_MISSING_VALUE)
    uniprot_reaction_map = (
        uniprot_df.set_index("uniprot_id")
        if "uniprot_id" in uniprot_df.columns
        else pd.DataFrame()
    )
    if "reaction_ec_numbers" in uniprot_reaction_map.columns:
        uniprot_reactions = resolved_for_combined.reindex(
            combined_df.index, fill_value=UNIPROT_MISSING_VALUE
        ).map(uniprot_reaction_map["reaction_ec_numbers"])
        uniprot_reactions = uniprot_reactions.fillna(UNIPROT_MISSING_VALUE)
    else:
        uniprot_reactions = pd.Series(
            UNIPROT_MISSING_VALUE,
            index=combined_df.index,
            dtype=object,
        )

    combined_df["reaction_ec_numbers"] = [
        normalize_reaction_ec_numbers([u, c])
        for u, c in zip(uniprot_reactions, chembl_reactions, strict=False)
    ]

    from tempfile import NamedTemporaryFile

    with NamedTemporaryFile(
        "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
    ) as tmp:
        combined_df.to_csv(
            tmp, index=False, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding
        )

        iuphar_input = Path(tmp.name)

    iuphar_args = argparse.Namespace(
        input_csv=iuphar_input, output_csv=output_csv, final_out=output_csv
    )
    orig_target = cfg.target.iuphar.target_csv
    orig_family = cfg.target.iuphar.family_csv
    cfg.target.iuphar.target_csv = cfg.target.all.target_csv
    cfg.target.iuphar.family_csv = cfg.target.all.family_csv
    try:
        if run_iuphar(cfg, iuphar_args) != 0:
            raise RuntimeError("IUPHAR retrieval failed")
    finally:
        cfg.target.iuphar.target_csv = orig_target
        cfg.target.iuphar.family_csv = orig_family
        iuphar_input.unlink(missing_ok=True)

    if not output_csv.exists():
        logger.warning(
            "missing_iuphar_output_file",
            path=str(output_csv),
        )
        _ensure_parent_directory(output_csv, cfg=cfg)
        empty_iuphar = pd.DataFrame({"uniprot_id": pd.Series(dtype=object)})
        write_csv_deterministic(
            empty_iuphar,
            output_csv,
            key_cols=["uniprot_id"],
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )

    try:
        iuphar_df = pd.read_csv(
            output_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
    except pd.errors.EmptyDataError:
        logger.warning("empty_iuphar_output", path=str(output_csv))
        iuphar_df = pd.DataFrame({"uniprot_id": pd.Series(dtype=object)})
    else:
        if "uniprot_id" not in iuphar_df.columns:
            logger.warning(
                "missing_iuphar_uniprot_column",
                path=str(output_csv),
                columns=list(iuphar_df.columns),
            )
            placeholder = pd.Series(
                UNIPROT_MISSING_VALUE, index=iuphar_df.index, dtype=object
            )
            iuphar_df = iuphar_df.copy()
            iuphar_df.insert(0, "uniprot_id", placeholder)
    combined_df, iuphar_df = _prepare_iuphar_merge_frames(combined_df, iuphar_df)
    logger.info("fetch_iuphar_done", rows=len(iuphar_df))
    return combined_df, iuphar_df


def merge_results(
    combined_df: pd.DataFrame,
    iuphar_df: pd.DataFrame,
    cfg: Config,
    *,
    classifier: ii.IUPHARClassifier | None = None,
) -> pd.DataFrame:
    """Merge ChEMBL, UniProt and IUPHAR data into a single table.

    Parameters
    ----------
    combined_df : pandas.DataFrame
        DataFrame containing ChEMBL and UniProt information produced by
        :func:`fetch_iuphar`.
    iuphar_df : pandas.DataFrame
        IUPHAR classification results aligned with ``combined_df`` on
        ``uniprot_id``.
    cfg : Config
        Application configuration providing classifier settings.
    classifier : library.integration.iuphar_library.IUPHARClassifier, optional
        Pre-initialised classifier. When ``None`` a classifier is created from
        ``cfg``.

    Returns
    -------
    pandas.DataFrame
        Merged and post-processed target data ready for validation and export.
    """

    logger.info("merge_results_start")
    chembl_keys: set[str] = set()
    if "uniprot_id" in combined_df.columns:
        chembl_keys = {
            str(value).strip()
            for value in combined_df["uniprot_id"].dropna()
            if str(value).strip()
        }
    iuphar_keys: set[str] = set()
    if "uniprot_id" in iuphar_df.columns:
        iuphar_keys = {
            str(value).strip()
            for value in iuphar_df["uniprot_id"].dropna()
            if str(value).strip()
        }
    logger.info(
        "merge_results_key_coverage",
        chembl=len(chembl_keys),
        iuphar=len(iuphar_keys),
        matched=len(chembl_keys & iuphar_keys),
    )
    merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
    if classifier is None:
        classifier = pc.classifier_from_config(cfg)
    merged = pc.append_protein_class_predictions(merged, classifier)
    processed = tp.postprocess_targets(merged)
    final_df = tp.finalise_targets(processed)
    logger.info("merge_results_done", rows=len(final_df))
    return final_df


def validate_and_write(
    df: pd.DataFrame,
    output: Path,
    cfg: Config,
    *,
    raw_out: Path | None = None,
    id_cols: Sequence[str] | None = None,
    raw_format: str = "csv",
    reindex_raw: bool = True,
    postprocess_context: IsoformPostprocessContext | None = None,
    table_name: str | None = None,
    date_tag: str | None = None,
    emit_legacy_artifacts: bool = False,
) -> int:
    """Normalise, validate and export the target table.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame produced by :func:`merge_results`.
    output : pathlib.Path
        Desired destination path used as a naming hint for the standard
        artefacts persisted via :func:`library.io.save_standard_outputs`.
    cfg : Config
        Application configuration providing schema definitions and IO settings.
    table_name : str, optional
        Explicit table identifier overriding the value inferred from
        ``output`` when generating filenames for the canonical artefacts.
    date_tag : str, optional
        Eight-digit UTC date tag used for deterministic artefact naming.
    emit_legacy_artifacts : bool, optional
        When ``True`` the legacy side outputs (raw exports, metadata sidecars)
        are produced alongside the standard CSV bundle.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate validation errors or
        failures when generating table-quality reports.
    """

    output_path = Path(output)
    inferred_table_name, inferred_date_tag = _resolve_output_metadata(
        output_path,
        date_hint=date_tag,
        table_hint=table_name,
    )
    io_cfg = cfg.io
    if hasattr(io_cfg, "model_copy"):
        output_io_cfg = io_cfg.model_copy(update={"output_dir": output_path.parent})
    elif hasattr(io_cfg, "copy"):
        output_io_cfg = io_cfg.copy(update={"output_dir": output_path.parent})  # type: ignore[arg-type]
    else:  # pragma: no cover - defensive fallback
        output_io_cfg = replace(io_cfg, output_dir=output_path.parent)
    expected_dataset_path = (
        Path(output_io_cfg.output_dir)
        / f"output.{inferred_table_name}_{inferred_date_tag}.csv"
    )
    logger.info("validate_write_start", output=str(expected_dataset_path))

    key_columns = list(id_cols) if id_cols else ["target_chembl_id"]
    input_rows = len(df)

    if emit_legacy_artifacts and raw_out is not None:
        raw_format_value = (raw_format or "csv").lower()
        raw_frame = df.copy()
        raw_order: Sequence[str] | None = TARGETS_COLUMN_ORDER if reindex_raw else None
        if raw_order is not None:
            raw_frame = raw_frame.reindex(columns=raw_order)
        if raw_format_value == "parquet":
            try:
                raw_frame.to_parquet(raw_out, index=False)
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ValueError(
                    "Parquet export requires optional pyarrow or fastparquet"
                ) from exc
        else:
            io.write_csv(
                raw_frame,
                raw_out,
                cfg=cfg,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
                key_cols=key_columns or None,
                col_order=raw_order,
            )
    elif raw_out is not None:
        logger.info(
            "legacy_raw_export_skipped",
            path=str(raw_out),
            reason="emit_legacy_artifacts_disabled",
        )

    normalized = normalize_targets(df)
    normalized_rows = len(normalized)
    logger.info(
        "normalize_targets_counts",
        before=input_rows,
        after=normalized_rows,
    )
    normalized = add_pipeline_metadata(normalized)
    prepared, missing_required, missing_optional = _prepare_targets_for_schema(
        normalized
    )
    final_df = prepared
    logger.info("prepare_targets_rows", rows=len(final_df))

    drop_reasons: dict[str, set[Any]] = {
        "schema": set(),
        "fk": set(),
        "regex": set(),
        "not_null": set(),
    }

    def _categorize_failure(row: pd.Series) -> str:
        check_value = str(row.get("check", "")).lower()
        context_value = str(row.get("schema_context", "")).lower()
        failure_value = str(row.get("failure_case", "")).lower()
        if "fk" in check_value or "foreign" in check_value or "fk" in context_value:
            return "fk"
        if "regex" in check_value or "match" in check_value or "pattern" in check_value:
            return "regex"
        if (
            "notnull" in check_value
            or "not null" in check_value
            or "nullable" in check_value
            or "null" in failure_value
        ):
            return "not_null"
        return "schema"

    exit_code = 0
    if not missing_required:
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )
        logger.info("targets_schema_validate_start", rows=len(final_df))
        try:
            validation = validate_targets(final_df, return_result=True)
        except SchemaErrors as exc:
            failure_path = expected_dataset_path.with_name(
                f"{expected_dataset_path.stem}_failure_cases.csv"
            )
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path, cfg=cfg)
            logger.error(
                "validation_failed",
                failures=len(exc.failure_cases),
                path=str(failure_path),
            )
            failure_cases = exc.failure_cases.copy()
            for _, failure_row in failure_cases.iterrows():
                reason = _categorize_failure(failure_row)
                identifier = failure_row.get("index")
                if pd.isna(identifier):
                    identifier = (
                        failure_row.get("column"),
                        failure_row.get("failure_case"),
                    )
                drop_reasons[reason].add(identifier)
            logger.info(
                "targets_schema_validate_result",
                status="failed",
                rows=len(final_df),
                failures=len(failure_cases),
            )
            final_df = getattr(exc, "validated_data", final_df)
            exit_code = 1
        else:
            final_df = validation.data
            failure_cases = validation.failure_cases.copy()
            if not failure_cases.empty:
                failure_path = expected_dataset_path.with_name(
                    f"{expected_dataset_path.stem}_failure_cases.csv"
                )
                errors = SidecarErrors()
                for row in failure_cases.to_dict("records"):
                    errors.add_error(row)
                errors.save(failure_path, cfg=cfg)
                logger.error(
                    "validation_failed",
                    failures=len(failure_cases),
                    path=str(failure_path),
                )
                for _, failure_row in failure_cases.iterrows():
                    reason = _categorize_failure(failure_row)
                    identifier = failure_row.get("index")
                    if pd.isna(identifier):
                        identifier = (
                            failure_row.get("column"),
                            failure_row.get("failure_case"),
                        )
                    drop_reasons[reason].add(identifier)
                logger.info(
                    "targets_schema_validate_result",
                    status="failed",
                    rows=len(final_df),
                    failures=len(failure_cases),
                )
                exit_code = 1
            else:
                logger.info(
                    "targets_schema_validate_result",
                    status="success",
                    rows=len(final_df),
                    failures=0,
                )
    else:

        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )

    total_dropped = sum(len(ids) for ids in drop_reasons.values())
    logger.info(
        "validation_drop_stats",
        total=total_dropped,
        schema=len(drop_reasons["schema"]),
        fk=len(drop_reasons["fk"]),
        regex=len(drop_reasons["regex"]),
        not_null=len(drop_reasons["not_null"]),
    )

    placeholders_before_fill = int(final_df.isna().sum().sum())
    logger.info("placeholder_fillna_pending", replacements=placeholders_before_fill)
    final_df = final_df.fillna("-")
    before_dedup = len(final_df)
    final_df = final_df.drop_duplicates()
    logger.info("deduplicated_rows", dropped=before_dedup - len(final_df))
    export_columns = [col for col in TARGETS_COLUMN_ORDER if col in final_df.columns]
    export_columns.extend(
        column for column in final_df.columns if column not in export_columns
    )
    final_df = final_df.loc[:, export_columns]
    if key_columns:
        sortable_keys = [col for col in key_columns if col in final_df.columns]
        if sortable_keys:
            final_df = final_df.sort_values(by=sortable_keys).reset_index(drop=True)
    quality_summary = generate_qc_report(
        final_df,
        table_name=inferred_table_name,
    )
    correlation_matrix = generate_correlation_report(
        final_df,
        table_name=inferred_table_name,
    )
    artifacts = io.save_standard_outputs(
        final_df,
        correlation_matrix,
        quality_summary,
        table_name=inferred_table_name,
        date_tag=inferred_date_tag,
        cfg=output_io_cfg,
    )
    logger.info(
        "standard_outputs_written",
        dataset=str(artifacts.dataset),
        quality=str(artifacts.quality_report),
        correlation=str(artifacts.correlation_report),
    )
    final_csv_path = artifacts.dataset
    ambiguous_count = 0
    if "protein_class_pred_rule_id" in final_df.columns:
        ambiguous_count = int(
            final_df["protein_class_pred_rule_id"]
            .fillna("")
            .astype(str)
            .str.lower()
            .eq("ambiguous")
            .sum()
        )
    http_requests = postprocess_context.http_requests if postprocess_context else None
    if http_requests is None:
        chembl_cfg = getattr(getattr(cfg, "target", None), "chembl", None)
        chunk_value = (
            getattr(chembl_cfg, "chunk_size", 0) if chembl_cfg is not None else 0
        )
        if chunk_value:
            http_requests = int(math.ceil(max(input_rows, 0) / chunk_value))
    isoform_context = IsoformPostprocessContext(
        args=postprocess_context.args if postprocess_context else None,
        http_requests=http_requests,
    )
    pipeline_result: "PostprocessingPipelineResult" | None = None
    if final_csv_path is not None:
        if exit_code != 0:
            logger.warning(
                "postprocess_running_with_errors",
                exit_code=exit_code,
                path=str(final_csv_path),
            )
        pipeline_result = run_target_postprocess_if_requested(
            final_csv_path,
            cfg=cfg,
            args=isoform_context.args,
            context=isoform_context,
            ambiguous_classifications=ambiguous_count,
        )
    if final_df.empty:
        logger.info(
            "quality_report_skipped",
            reason="empty_dataframe",
            table=str(final_csv_path.with_suffix("")),
        )
    else:
        doc_quality_cfg = cfg.system.doc_quality
        table_quality = build_table_quality_hook(
            doc_quality_cfg,
            table_name=final_csv_path.with_suffix(""),
            destination=final_csv_path.parent,
        )
        try:
            table_quality(final_df)
        except Exception as exc:
            logger.exception(
                "quality_report_failed",
                error=str(exc),
                path=str(final_csv_path),
                exc=exc,
            )
            return 1

    metrics = pipeline_result.metrics if pipeline_result else None
    report_path = pipeline_result.report_path if pipeline_result else None
    pipeline_version_value = (
        metrics.pipeline_version
        if metrics and metrics.pipeline_version is not None
        else get_pipeline_version()
    )
    payload: dict[str, object] = {
        "rows": len(final_df),
        "output_postprocessed": str(final_csv_path),
        "pipeline_version": pipeline_version_value,
    }
    payload["quality_report"] = str(artifacts.quality_report)
    payload["correlation_report"] = str(artifacts.correlation_report)
    if metrics is not None:
        summary = metrics.summary()
        if summary.get("rows") is not None:
            payload["postprocess_rows"] = summary["rows"]
        if summary.get("columns") is not None:
            payload["postprocess_columns"] = summary["columns"]
        if summary.get("duration_s") is not None:
            payload["postprocess_duration_s"] = summary["duration_s"]
        if summary.get("steps") is not None:
            payload["postprocess_steps"] = summary["steps"]
        if metrics.validation is not None:
            payload["postprocess_schema"] = metrics.validation.schema
    if report_path is not None:
        payload["postprocess_report"] = str(report_path)
    logger.info("validate_write_done", **payload)
    return exit_code


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run the full target acquisition pipeline.

    Parameters
    ----------
    cfg : Config
        Application configuration combining defaults for every individual
        pipeline.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero values indicate failures while fetching
        data, validating the merged dataset or writing the resulting CSV. The
        offending step is logged in the ``pipeline_step_failed`` event.
    """

    limit = cfg.target.all.limit
    if limit is not None and limit < 1:
        logger.error(
            "invalid_limit",
            section="target.all.limit",
            limit=limit,
        )
        return 1

    disable_gtop_cli = bool(getattr(args, "disable_gtop", False))
    original_enable_gtop = cfg.target.uniprot.enable_gtop
    if disable_gtop_cli:
        cfg.target.uniprot.enable_gtop = False
    try:
        final_candidate = getattr(args, "final_out", None)
        if final_candidate in (None, argparse.SUPPRESS):
            final_output = Path(
                io.default_output_path(
                    args.input_csv,
                    cfg.io,
                    date=getattr(args, "date", None),
                )
            )
        else:
            final_output = Path(final_candidate)

        emit_legacy = bool(getattr(args, "emit_legacy_artifacts", False))
        raw_candidate = getattr(args, "raw_out", None)
        if not emit_legacy or raw_candidate in (None, argparse.SUPPRESS):
            raw_output: Path | None = None
        else:
            raw_output = Path(raw_candidate)
        chembl_out = cfg.target.all.chembl_out or final_output.with_name(
            final_output.stem + "_chembl.csv"
        )
        uniprot_out = cfg.target.all.uniprot_out or final_output.with_name(
            final_output.stem + "_uniprot.csv"
        )
        iuphar_out = cfg.target.all.iuphar_out or final_output.with_name(
            final_output.stem + "_iuphar.csv"
        )

        raw_format = str(getattr(args, "raw_format", "csv") or "csv").lower()
        reindex_raw = not bool(getattr(args, "no_reindex_raw", False))
        id_cols_value = getattr(args, "id_cols", None)
        if id_cols_value in (None, argparse.SUPPRESS):
            key_columns = ["target_chembl_id"]
        elif isinstance(id_cols_value, str):
            key_columns = [id_cols_value]
        else:
            key_columns = list(id_cols_value) or ["target_chembl_id"]

        table_hint = getattr(args, "_auto_output_stem", None)
        raw_date_hint = getattr(args, "date", None)
        date_hint = raw_date_hint.strip() if isinstance(raw_date_hint, str) else None
        resolved_table_name, resolved_date_tag = _resolve_output_metadata(
            final_output,
            date_hint=date_hint,
            table_hint=table_hint,
        )

        chembl_df = fetch_chembl(
            cfg,
            args.input_csv,
            chembl_out,
            limit=limit,
            chunk_size=cfg.target.all.chunk_size,
            offset=cfg.target.all.offset,
            id_cols=key_columns,
        )
        uniprot_df = fetch_uniprot(cfg, chembl_df, uniprot_out)
        combined_df, iuphar_df = fetch_iuphar(cfg, chembl_df, uniprot_df, iuphar_out)
        merged = merge_results(combined_df, iuphar_df, cfg)
        exit_code = validate_and_write(
            merged,
            final_output,
            cfg,
            raw_out=raw_output,
            id_cols=key_columns,
            raw_format=raw_format,
            reindex_raw=reindex_raw,
            postprocess_context=IsoformPostprocessContext(args=args),
            table_name=resolved_table_name,
            date_tag=resolved_date_tag,
            emit_legacy_artifacts=emit_legacy,
        )
        return exit_code
    except (FileNotFoundError, ValueError, OSError, RuntimeError) as exc:
        logger.error(
            "pipeline_step_failed",
            error=str(exc),
            exc_info=exc,
            step="all",
            input=str(args.input_csv),
            output=str(final_output),
        )
        return 1
    finally:
        cfg.target.uniprot.enable_gtop = original_enable_gtop


def _update_target_config_from_options(
    cfg: Config, options: TargetPipelineOptions
) -> None:
    """Apply limit/offset overrides from ``options`` to ``cfg``."""

    target_cfg = cfg.sources.chembl.pipelines.target

    def _apply(section_name: str) -> None:
        section = getattr(target_cfg, section_name)
        updates: dict[str, object] = {"offset": options.offset}
        if options.limit is not None:
            updates["limit"] = options.limit
        setattr(target_cfg, section_name, section.model_copy(update=updates))

    command = options.command
    if command == "all":
        for section_name in ("all", "chembl", "uniprot", "iuphar"):
            _apply(section_name)
    elif command in {"chembl", "uniprot", "iuphar"}:
        _apply(command)
    else:  # pragma: no cover - defensive guard
        msg = f"unsupported target pipeline command: {command}"
        raise ValueError(msg)


def run_target_service(
    config: Config, options: TargetPipelineOptions
) -> PipelineRunResult:
    """Execute the target pipeline using typed ``options``."""

    output_path = Path(options.output_csv)
    if options.skip_existing and output_path.exists() and not options.force:
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=False,
            reason="skip_existing",
            written=False,
        )

    cfg = config.model_copy(deep=True)
    _update_target_config_from_options(cfg, options)

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        final_out=output_path,
        output_csv=output_path,
        raw_out=options.raw_output,
        raw_format=options.raw_format,
        no_reindex_raw=options.no_reindex_raw,
        id_cols=options.id_columns,
        limit=options.limit,
        offset=options.offset,
        command=options.command,
        skip_existing=False,
        force=options.force,
    )

    command = options.command
    if command == "chembl":
        runner = run_chembl
    elif command == "uniprot":
        runner = run_uniprot
    elif command == "iuphar":
        runner = run_iuphar
    elif command == "all":
        runner = run_all
    else:  # pragma: no cover - defensive guard
        msg = f"unsupported target pipeline command: {command}"
        raise ValueError(msg)

    exit_code = int(runner(cfg, args))
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the selected target pipeline applying CLI policies."""

    alias_token = getattr(args, "_command_keyboard_alias", None)
    if isinstance(alias_token, str):
        logger.info(
            "command_alias_resolved",
            alias=alias_token,
            command=args.command,
        )

    parameter_entries = _resolve_target_parameters(args.command, cfg, args)
    for entry in parameter_entries:
        logger.info(
            "cli_parameter",
            command=args.command,
            scope=entry.scope,
            param=entry.name,
            value=entry.value,
            source=entry.source,
        )

    final_candidate = getattr(args, "final_out", None)
    if final_candidate in (None, argparse.SUPPRESS):
        final_output = Path(
            io.default_output_path(
                args.input_csv,
                cfg.io,
                date=getattr(args, "date", None),
            )
        )
        args.final_out = final_output
    else:
        final_output = Path(final_candidate)
        if not isinstance(final_candidate, Path):
            args.final_out = final_output
        else:
            args.final_out = final_candidate
    if args.skip_existing and final_output.exists() and not args.force:
        logger.info("pipeline_skip_existing", output=str(final_output))
        return 0
    func = getattr(args, "func", None)
    if func is None:
        logger.error("missing_subcommand_handler", command=getattr(args, "command", ""))
        return 1
    result = func(cfg, args)
    return int(result)


class TargetPipelineCLI(PipelineCLIBase):
    """CLI adapter orchestrating the target data pipelines."""

    def build_parser(self) -> tuple[argparse.ArgumentParser, LoggerConfig]:
        return _build_parser_impl()

    def prepare_arguments(
        self,
        parser: argparse.ArgumentParser,
        args: argparse.Namespace,
        argv: Sequence[str] = None,
    ) -> argparse.Namespace:
        del argv
        command_alias = COMMAND_ALIAS_TO_CANONICAL.get(getattr(args, "command", ""))
        if command_alias is not None:
            args._command_keyboard_alias = args.command
            args.command = command_alias
        prepare_io_paths(
            args,
            input_default=DEFAULT_INPUT_NAME,
            output_stem=DEFAULT_OUTPUT_STEM,
        )
        return args

    def get_program_name(self) -> str:
        return Path(__file__).with_suffix("").name

    def get_logger(self) -> Logger:
        return logger

    def execute(
        self,
        args: argparse.Namespace,
        parser: argparse.ArgumentParser,
        logging_ctx: CLILoggingContext,
    ) -> int:
        configure_logger(logging_ctx.log_cfg)
        console_stream = logging_ctx.console_stream
        exit_code = 0
        try:
            limit_value = getattr(args, "limit", None)
            subparser_map = getattr(parser, "subparsers_map", {})
            subparser = subparser_map.get(getattr(args, "command", None), parser)
            if limit_value is not None and limit_value < 1:
                subparser.error("--limit must be a positive integer")
            offset_value = getattr(args, "offset", 0)
            if offset_value < 0:
                subparser.error("--offset must be zero or a positive integer")
            mapping: dict[str, str] = {}
            args_namespace: argparse.Namespace = args
            if args.command == "uniprot":
                mapping = {
                    "column": "target.uniprot.column",
                    "chunk_size": "target.uniprot.chunk_size",
                    "timeout": "target.uniprot.timeout",
                    "limit": "target.uniprot.limit",
                    "offset": "target.uniprot.offset",
                    "data_dir": "target.uniprot.data_dir",
                }
            elif args.command == "chembl":
                mapping = {
                    "column": "target.chembl.column",
                    "chunk_size": "target.chembl.chunk_size",
                    "timeout": "target.chembl.timeout",
                    "limit": "target.chembl.limit",
                    "offset": "target.chembl.offset",
                }
            elif args.command == "iuphar":
                mapping = {
                    "column": "target.iuphar.column",
                    "chunk_size": "target.iuphar.chunk_size",
                    "timeout": "target.iuphar.timeout",
                    "limit": "target.iuphar.limit",
                    "offset": "target.iuphar.offset",
                    "target_csv": "target.iuphar.target_csv",
                    "family_csv": "target.iuphar.family_csv",
                }
            elif args.command == "all":
                mapping = {
                    "column": "target.all.column",
                    "chunk_size": "target.all.chunk_size",
                    "timeout": "target.all.timeout",
                    "limit": "target.all.limit",
                    "offset": "target.all.offset",
                    "data_dir": "target.all.data_dir",
                    "target_csv": "target.all.target_csv",
                    "family_csv": "target.all.family_csv",
                    "uniprot_column": "target.all.uniprot_column",
                    "chembl_out": "target.all.chembl_out",
                    "uniprot_out": "target.all.uniprot_out",
                    "iuphar_out": "target.all.iuphar_out",
                    "chembl_column": "target.chembl.column",
                    "chembl_chunk_size": "target.chembl.chunk_size",
                    "chembl_timeout": "target.chembl.timeout",
                    "chembl_limit": "target.chembl.limit",
                    "chembl_offset": "target.chembl.offset",
                    "uniprot_chunk_size": "target.uniprot.chunk_size",
                    "uniprot_timeout": "target.uniprot.timeout",
                    "uniprot_limit": "target.uniprot.limit",
                    "uniprot_offset": "target.uniprot.offset",
                    "iuphar_column": "target.iuphar.column",
                    "iuphar_chunk_size": "target.iuphar.chunk_size",
                    "iuphar_timeout": "target.iuphar.timeout",
                    "iuphar_limit": "target.iuphar.limit",
                    "iuphar_offset": "target.iuphar.offset",
                }
                args_dict = vars(args).copy()
                final_value = args_dict.get("final_out")
                if final_value in (None, argparse.SUPPRESS):
                    temp_namespace = argparse.Namespace(**args_dict)
                    cli.prepare_io_paths(
                        temp_namespace,
                        input_default=DEFAULT_INPUT_NAME,
                        output_stem=DEFAULT_OUTPUT_STEM,
                    )
                    args_dict["final_out"] = temp_namespace.final_out
                    args_dict["output_csv"] = temp_namespace.output_csv
                    args_dict["raw_out"] = temp_namespace.raw_out
                    args_dict["date"] = getattr(temp_namespace, "date", None)
                else:
                    args_dict["final_out"] = Path(final_value)
                args_namespace = argparse.Namespace(**args_dict)
            if mapping:
                exit_code = run_cli_command(
                    args=args_namespace,
                    parser=subparser,
                    base_parser=parser,
                    log_cfg=logging_ctx.log_cfg,
                    mapping=mapping,
                    run=run,
                    logger=logger,
                )
        except PipelineError as exc:
            exit_code = 2
            logger.error("pipeline_error", error=str(exc), exc_info=exc)
            print(f"[ERROR] {exc}", file=console_stream, flush=True)
        return exit_code


_CLI = TargetPipelineCLI()


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Expose parser creation for backwards compatibility."""

    return _CLI.build_parser()


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate to :class:`TargetPipelineCLI` for backwards compatibility."""

    return _CLI.main(argv)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
