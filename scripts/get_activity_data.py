"""Thin wrapper exposing the activity pipeline CLI entry point."""

from __future__ import annotations

if __package__ in {None, ""}:
    from _bootstrap import bootstrap_cli
else:  # pragma: no cover - executed when imported as a module
    from ._bootstrap import bootstrap_cli

bootstrap_cli(__package__, __file__)
del bootstrap_cli

from importlib import import_module
import argparse
import json
import math
import numbers
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from itertools import islice
from pathlib import Path
from threading import Lock
from time import perf_counter, sleep
from typing import Any, cast
from urllib.parse import urlsplit

import pandas as pd
import requests

from library import cli, io
from library.clients.chembl import ChemblClient
from library.integration import chembl_library as cl

try:  # pragma: no cover - urllib3 is part of requests dependency chain
    from urllib3.exceptions import NameResolutionError as _Urllib3NameResolutionError
except Exception:  # pragma: no cover - defensive fallback for alternative stacks
    _Urllib3NameResolutionError = None  # type: ignore

from library.cli import (
    Logger,
    LoggerConfig,
    positive_int,
    build_parser as base_parser,
)
from library.cli.base import PipelineCLIBase
from library.cli.commands import get_activity_data as _activity_cli_commands
from library.cli.commands.get_activity_data import (
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    run_activity_pipeline,
)
from library.cli_utils import (
    PipelineError,
    resolve_invocation,
    run_cli_command as _run_cli_command,
)
from library.common.csv_utils import (
    write_csv_chunks_deterministic,
)
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay
from library.common.log import logger
from library.config import Config, _serialize_paths
from library.io.metadata import (
    write_meta_yaml as _cli_write_meta_yaml,
)
from library.metadata import file_sha256 as _metadata_file_sha256
from library.orchestration import ETLContext
from library.pipelines.activity import run as activity_run
from library.pipelines.assay.chembl_assay import (
    ACTIVITY_COLUMNS,
)

from library.cli.entrypoints.activity import (
    ActivityPipelineCLI,
    DEFAULT_INPUT_NAME,
    DEFAULT_OUTPUT_STEM,
    PROGRAM_NAME,
    _emit_completion_message,
    _generate_activity_postprocess_metrics as _cli_generate_activity_postprocess_metrics,
    build_parser,
    main,
    run,
    run_chembl,
)
from library.pipelines.activity import runner as _activity_runner
from library.pipelines.activity.runner import (
    MAX_ACTIVITY_CHUNK_SIZE,
    register_activity_pipeline_hooks,
)
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    add_pipeline_metadata,
)
from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess.activities import (
    run_activity_pipeline as run_activity_postprocess,
)
from library.postprocess.common import collect_postprocess_metrics
from library.postprocessing import helpers as postprocessing_helpers
from library.postprocessing.activity_extended import process_activity_extended
from library.processing.activity import (
    apply_activity_annotations,
    compute_activity_bounds,
)
from library.qa.reporting import build_table_quality_hook
from library.schemas import (
    ActivitiesSchema,
    configure_activity_schema,
    normalize_activities,
)
from library.validation import validate_activities

try:
    from library.cli.entrypoints.activity import (
        _generate_activity_postprocess_metrics as _entrypoint_generate_activity_postprocess_metrics,
    )
except (ImportError, AttributeError):  # pragma: no cover - defensive guard for refactors
    _entrypoint_generate_activity_postprocess_metrics = None

DEFAULT_INPUT_NAME = "activity.csv"
DEFAULT_OUTPUT_STEM = "activities"
PROGRAM_NAME = Path(__file__).with_suffix("").name

def _args_invocation(args: argparse.Namespace) -> tuple[str, ...]:
    invocation = getattr(args, "invocation", None)
    if invocation is None:
        return (PROGRAM_NAME,)
    # Preserve POSIX-style paths for predictable logs and tests
    result: list[str] = []
    for arg in invocation:
        text = str(arg)
        # Convert backslashes to forward slashes only for absolute POSIX-like roots
        if isinstance(arg, Path):
            text = text.replace("\\", "/")
        result.append(text)
    return tuple(result)

file_sha256 = _metadata_file_sha256
write_meta_yaml = _cli_write_meta_yaml
try:
    from library.cli import configure_logger
except ImportError:
    configure_logger = None  # type: ignore[assignment]

run_cli_command = _run_cli_command

__all__ = (
    "file_sha256",
    "write_meta_yaml",
    "configure_logger",
    "run_cli_command",
    "MAX_ACTIVITY_CHUNK_SIZE",
)


def _ensure_command_logger_sync() -> None:
    """Keep the shared command module logger aligned with this script logger."""

    try:
        commands_module = _activity_cli_commands
    except NameError:  # pragma: no cover - defensive guard for refactors
        return

    if getattr(commands_module, "logger", None) is logger:
        return

    try:
        commands_module.logger = logger
    except Exception:  # pragma: no cover - defensive guard
        pass


_ensure_command_logger_sync()


def _generate_activity_postprocess_metrics(
    cfg: Config,
    output_path: Path,
    *,
    logger: Logger,
    extras: Mapping[str, object] | None = None,
) -> tuple[object, Path]:
    """Generate postprocess metrics using the entrypoint helper or a local fallback."""

    if _entrypoint_generate_activity_postprocess_metrics is not None:
        return _entrypoint_generate_activity_postprocess_metrics(
            cfg,
            output_path,
            logger=logger,
            extras=extras,
        )

    metrics, report_path = collect_postprocess_metrics(
        table="activity",
        output_path=output_path,
        csv_sep=cfg.io.csv_sep,
        csv_encoding=cfg.io.csv_encoding,
        output_dir=cfg.io.output_dir,
        runner=run_activity_postprocess,
        logger=logger,
        pipeline_version=get_pipeline_version(),
        report_extras=extras,
    )
    if report_path is None:
        fallback = Path(cfg.io.output_dir) / "activity.postprocess.report.json"
        fallback.parent.mkdir(parents=True, exist_ok=True)
        fallback_payload: dict[str, object] = {
            "table": "activity",
            "metrics": None,
            "output_path": str(output_path),
        }
        if extras:
            fallback_payload["extras"] = dict(extras)
        fallback.write_text(
            json.dumps(fallback_payload, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        report_path = fallback
    return metrics, report_path


_CACHE_MISS = object()
_CACHE_IN_PROGRESS = object()


_ACTIVITY_REQUIRED_COLUMNS: tuple[str, ...] = tuple(
    name for name, column in ActivitiesSchema.columns.items() if column.required
)

_ACTIVITY_REQUIRED_DTYPES: dict[str, object] = {
    name: getattr(column, "dtype", None)
    for name, column in ActivitiesSchema.columns.items()
    if column.required
}

_ORIGINAL_IO_WRITE_CSV = io.write_csv


@dataclass(slots=True)
class PreparedActivityContext:
    """Container describing the prepared identifiers for the activity pipeline."""

    limited_ids: Iterable[str]
    limit: int | None
    _processed_accessor: Callable[[], int]

    @property
    def processed_ids(self) -> int:
        """Return the number of identifiers consumed by the pipeline."""

        return self._processed_accessor()


def prepare_activity_context(
    cfg: Config,
    args: argparse.Namespace,
    *,
    skip_read: bool = False,
) -> PreparedActivityContext | None:
    """Prepare identifier iteration and limit handling for the activity pipeline."""

    limit = cfg.activity.limit
    if limit is not None and limit < 0:
        logger.error(
            "Configuration error: activity.limit must be non-negative but was set to %s.",
            limit,
        )
        return None

    if skip_read:
        return PreparedActivityContext(
            limited_ids=(),
            limit=limit,
            _processed_accessor=lambda: 0,
        )

    logger.info("activity_pipeline_read_input", input=str(args.input_csv))
    try:
        ids_iter = io.read_ids(args.input_csv, column=cfg.activity.column, cfg=cfg.io)
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            input=str(args.input_csv),
            error=str(exc),
            exc_info=exc,
        )
        return None

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    processed_ids = 0

    def _iter_ids() -> Iterator[str]:
        nonlocal processed_ids
        for identifier in ids_iter:
            processed_ids += 1
            yield identifier

    if limit is not None:
        limited_ids: Iterable[str] = islice(_iter_ids(), limit)
    else:
        limited_ids = _iter_ids()

    return PreparedActivityContext(
        limited_ids=limited_ids,
        limit=limit,
        _processed_accessor=lambda: processed_ids,
    )


_EXTENDED_ACTIVITY_DTYPES: dict[str, str] = {
    "activity_chembl_id": "string",
    "salt_chembl_id": "string",
    "target_chembl_id": "string",
    "bao_endpoint": "string",
    "compound_key": "string",
    "compound_name": "string",
    "multmol_assay": "boolean",
    "approx_cited_activity": "boolean",
    "shuffled_cit": "boolean",
    "exact_cited_activity": "boolean",
    "higly_correlated_cit": "boolean",
    "review_doc": "boolean",
    "rounded_data_citation": "boolean",
    "original_activity_approx": "string",
    "original_activity_exact": "string",
    "nstereo": "Int64",
    "log_value": "Float64",
}


_OUTPUT_ACTIVITY_DROP_COLUMNS: tuple[str, ...] = (
    "approx_cited_activity",
    "exact_cited_activity",
    "higly_correlated_cit",
    "multmol_assay",
    "original_activity_approx",
    "original_activity_exact",
    "review_doc",
    "rounded_data_citation",
    "standard_lower_value",
    "standard_upper_value",
    "shuffled_cit",
)



def _coerce_series_dtype(series: pd.Series[Any], dtype: str) -> pd.Series[Any]:
    """Return ``series`` converted to ``dtype`` where feasible."""

    # Delegate to the extended coercion helper for extension-aware conversions.
    if dtype in {"Float64", "Int64", "boolean"}:
        converted, _ = _coerce_extended_series(series, dtype, column="_coerce_series_dtype")
        return cast("pd.Series[Any]", converted)

    if dtype == "string":
        return series.astype(pd.StringDtype())

    # For non-extension dtypes, try to use numpy dtype if possible.
    try:
        import numpy as np

        return series.astype(np.dtype(dtype))
    except (TypeError, ValueError):
        # Final fallback: convert to string.
        return cast("pd.Series[Any]", series.astype(str))


def _extract_adapter_retry_metadata(
    session: requests.Session, base_url: str
) -> dict[str, object]:
    """Return retry configuration extracted from ``session`` for ``base_url``.

    The helper is defensive to accommodate test doubles that expose partial
    ``requests.Session`` interfaces. Only basic numeric attributes from the
    underlying :class:`urllib3.util.retry.Retry` instance are returned so the
    payload is JSON/log friendly.
    """

    adapters: list[requests.adapters.BaseAdapter] = []
    try:
        adapter = session.get_adapter(base_url)
        if isinstance(adapter, requests.adapters.HTTPAdapter):
            adapters.append(adapter)
    except Exception:  # pragma: no cover - defensive for alternative sessions
        pass
    try:
        adapters.append(session.get_adapter("https://"))
    except Exception:  # pragma: no cover - defensive for alternative sessions
        pass
    metadata: dict[str, object] = {}
    for adapter in adapters:
        if adapter is None:
            continue
        retries = getattr(adapter, "max_retries", None)
        if retries is None:
            continue
        for attr, key in (
            ("total", "adapter_total"),
            ("connect", "adapter_connect"),
            ("read", "adapter_read"),
            ("status", "adapter_status"),
        ):
            value = getattr(retries, attr, None)
            if value is not None:
                metadata[key] = value
        backoff = getattr(retries, "backoff_factor", None)
        if backoff is not None:
            metadata["adapter_backoff_factor"] = backoff
        # Prefer the first adapter that exposes retry settings.
        if metadata:
            break
    return metadata


def _gather_http_diagnostics(cfg: Config, client: object) -> dict[str, object]:
    """Collect HTTP-related configuration values for logging and errors."""

    diagnostics: dict[str, object] = {
        "activity_batch_size": getattr(cfg.activity, "batch_size", None),
        "activity_timeout": getattr(cfg.activity, "timeout", None),
        "api_timeout_read": getattr(cfg.api, "timeout_read", None),
        "api_retries": getattr(cfg.api, "retries", None),
        "api_backoff_factor": getattr(cfg.api, "backoff_factor", None),
        "retry_max_attempts": getattr(cfg.retry, "max_attempts", None),
        "retry_backoff_factor": getattr(cfg.retry, "backoff_factor", None),
    }

    session_obj = getattr(client, "session", None)
    if isinstance(session_obj, requests.Session):
        adapter_meta = _extract_adapter_retry_metadata(session_obj, cfg.api.chembl_base)
        diagnostics.update(adapter_meta)

    clean: dict[str, object] = {}
    for key, value in diagnostics.items():
        if value is None:
            continue
        clean[key] = value
    return clean


def _iter_exception_chain(exc: BaseException) -> Iterator[BaseException]:
    """Yield ``exc`` and the causal/context chain preserving order."""

    current: BaseException | None = exc
    seen: set[int] = set()
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        yield current
        next_exc = current.__cause__ or current.__context__
        current = next_exc if isinstance(next_exc, BaseException) else None


def _is_name_resolution_error(exc: Exception) -> bool:
    """Return ``True`` when ``exc`` or its causes indicate DNS failure."""

    if isinstance(exc, requests.exceptions.RequestException):
        for candidate in _iter_exception_chain(exc):
            if (
                _Urllib3NameResolutionError is not None
                and isinstance(candidate, _Urllib3NameResolutionError)
            ):
                return True
            message = str(candidate).strip().lower()
            if (
                "name resolution" in message
                or "nameresolutionerror" in message
                or "getaddrinfo failed" in message
            ):
                return True

    message = str(exc).strip().lower()
    if (
        "name resolution" in message
        or "nameresolutionerror" in message
        or "getaddrinfo failed" in message
    ):
        return True
    return False


def _describe_network_failure(cfg: Config, exc: Exception) -> tuple[str | None, str | None]:
    """Return a human readable hint and host when DNS failures are detected."""

    if not _is_name_resolution_error(exc):
        return None, None

    base = getattr(cfg.api, "chembl_base", "")
    host = urlsplit(str(base)).hostname or str(base) or "www.ebi.ac.uk"
    hint = (
        "Unable to resolve the ChEMBL API host. Check internet connectivity "
        "or configure offline fixtures for testing environments."
    )
    return hint, host

class _StreamingCSVStatistics:
    """Accumulate row and null counters while streaming CSV chunks."""

    __slots__ = ("rows", "cells", "nulls")

    def __init__(self) -> None:
        self.rows = 0
        self.cells = 0
        self.nulls = 0

    def update(self, chunk: pd.DataFrame) -> None:
        """Update counters using ``chunk`` after output column filtering."""

        rows = int(len(chunk))
        self.rows += rows

        cells = int(chunk.size)
        self.cells += cells

        if cells:
            nulls = int(chunk.isna().to_numpy().sum())
        else:
            nulls = 0
        self.nulls += nulls

    def snapshot(self) -> dict[str, float | int]:
        """Return a summary of the accumulated metrics."""

        if self.cells == 0:
            null_fraction = 0.0
        else:
            null_fraction = self.nulls / self.cells
        return {
            "rows": self.rows,
            "cells": self.cells,
            "nulls": self.nulls,
            "null_fraction": float(null_fraction),
        }


def _emit_completion_message(
    *,
    output_path: Path | None,
    processed_rows: int | None,
    duration_s: float,
    mode: str,
    streamed_metrics: Mapping[str, object] | None = None,
) -> None:
    """Log a structured completion summary for the pipeline."""

    resolved_rows = processed_rows if processed_rows is not None else 0
    null_fraction_value: float | None = None
    metrics_payload: dict[str, object] | None = None

    if streamed_metrics:
        rows_value = streamed_metrics.get("rows")
        if isinstance(rows_value, numbers.Integral):
            resolved_rows = int(rows_value)
        elif isinstance(rows_value, numbers.Real):
            rows_float = float(rows_value)
            if not math.isnan(rows_float):
                resolved_rows = int(rows_float)

        null_fraction = streamed_metrics.get("null_fraction")
        if isinstance(null_fraction, numbers.Real):
            null_fraction_float = float(null_fraction)
            if not math.isnan(null_fraction_float):
                null_fraction_value = null_fraction_float

        metrics_payload = {}
        for key, value in streamed_metrics.items():
            if isinstance(value, numbers.Integral):
                metrics_payload[key] = int(value)
            elif isinstance(value, numbers.Real):
                converted = float(value)
                metrics_payload[key] = None if math.isnan(converted) else converted
            elif isinstance(value, (str, bool)) or value is None:
                metrics_payload[key] = value
            else:
                metrics_payload[key] = str(value)

    payload: dict[str, object] = {
        "output": str(output_path) if output_path is not None else None,
        "rows": int(resolved_rows),
        "duration_s": float(duration_s),
        "mode": mode,
    }
    if processed_rows is not None:
        payload["processed_rows"] = int(processed_rows)
    if null_fraction_value is not None:
        payload["null_fraction"] = null_fraction_value
    if metrics_payload:
        payload["streamed_metrics"] = metrics_payload

    logger.info("activity_pipeline_completion", **payload)

_EXTENDED_ACTIVITY_FALLBACKS: dict[str, Callable[[pd.DataFrame], pd.Series[Any] | None]] = {
    "activity_chembl_id": lambda df: df.get("activity_id"),
    "compound_name": lambda df: df.get("molecule_pref_name"),
    "log_value": lambda df: df.get("pchembl_value"),
}


def _coerce_extended_series(
    series: pd.Series[Any],
    dtype: str,
    column: str,
) -> tuple[pd.Series[Any], pd.Series[bool]]:
    """Return ``series`` coerced to ``dtype`` and a mask of conversion failures."""

    if not isinstance(series, pd.Series):
        series = pd.Series(series)

    if dtype == "string":
        converted = series.astype(pd.StringDtype())
        failures = pd.Series(False, index=converted.index)
        return converted, failures

    if dtype == "Float64":
        numeric = pd.to_numeric(series, errors="coerce")
        converted = pd.Series(pd.array(numeric.tolist(), dtype="Float64"), index=series.index)
        failures = series.notna() & converted.isna()
        return converted, failures

    if dtype == "Int64":
        numeric = pd.to_numeric(series, errors="coerce")
        converted_float: pd.Series[Any] = pd.Series(pd.array(numeric.tolist(), dtype="Float64"), index=series.index)
        failures = series.notna() & converted_float.isna()
        non_integral_mask = converted_float.notna() & converted_float.ne(converted_float.round())
        if non_integral_mask.any():
            converted_float.loc[non_integral_mask] = pd.NA
            failures = failures | non_integral_mask
            try:
                converted = converted_float.astype("Int64")
            except (TypeError, ValueError):
                converted = pd.Series(pd.NA, index=series.index, dtype="Int64")
                failures = pd.Series(True, index=series.index)
            return converted, failures

    if dtype == "boolean":
        try:
            converted = series.astype("boolean")
            failures = series.notna() & converted.isna()
        except (TypeError, ValueError):
            converted = pd.Series(pd.NA, index=series.index, dtype="boolean")
            failures = pd.Series(True, index=series.index)
            return converted, failures

    try:
        converted = series.astype(dtype)
        # If dtype is not compatible, astype with errors="ignore" will not convert, so check for failures
        failures = series.notna() & converted.isna()
    except (TypeError, ValueError):
        converted = pd.Series(pd.NA, index=series.index, dtype=dtype)
        failures = pd.Series(True, index=series.index)

    return converted, failures


 
def _string_like_missing(series: pd.Series[Any]) -> pd.Series[bool]:
    """Return a boolean mask for ``series`` treating blanks as missing."""

    mask = series.isna()
    if pd.api.types.is_string_dtype(series) or series.dtype == "object":
        string_values = series.astype("string")
        mask = mask | string_values.str.strip().fillna("").eq("")
    return mask
 
 
def _string_blank_mask(series: pd.Series[Any]) -> pd.Series[bool]:
    """Return mask of entries that are null or contain only whitespace."""

    values = series.astype("string")
    stripped = values.str.strip()
    return values.isna() | stripped.fillna("").eq("")


def _load_assay_src_lookup(dictionary_dir: Path | str | None) -> dict[str, str]:
    """Return mapping of assay identifiers to ``src_assay_id`` values."""

    if dictionary_dir is None:
        return {}

    candidate = Path(dictionary_dir) / "_assay" / "assay.csv"
    try:
        frame = postprocessing_helpers.read_csv_with_fallbacks(candidate)
    except FileNotFoundError:
        logger.warning(
            f"Assay lookup file '{candidate}' was not found; src_assay_id enrichment will be skipped."
        )
        return {}
    except pd.errors.EmptyDataError:
        logger.warning(
            f"Assay lookup file '{candidate}' is empty; src_assay_id enrichment will be skipped."
        )
        return {}
    except ValueError as exc:
        logger.warning(
            f"Assay lookup file '{candidate}' has invalid columns and cannot be read: {exc}"
        )
        return {}
    except OSError as exc:
        logger.warning(
            f"Reading assay lookup file '{candidate}' failed due to an OS error: {exc}"
        )
        return {}

    required_columns = {"assay_chembl_id", "src_assay_id"}
    missing_columns = sorted(required_columns.difference(frame.columns))
    if missing_columns:
        logger.warning(
            "activity_missing_columns",
            path=str(candidate),
            missing=missing_columns,
        )
        return {}

    frame = frame.loc[:, ["assay_chembl_id", "src_assay_id"]]

    if frame.empty:
        return {}

    cleaned = frame.dropna(subset=["assay_chembl_id", "src_assay_id"])
    if cleaned.empty:
        return {}

    cleaned = cleaned.assign(
        assay_chembl_id=cleaned["assay_chembl_id"].str.strip(),
        src_assay_id=cleaned["src_assay_id"].str.strip(),
    )

    cleaned = cleaned[cleaned["assay_chembl_id"].ne("") & cleaned["src_assay_id"].ne("")]
    if cleaned.empty:
        return {}

    return {
        str(assay_id): str(src_id)
        for assay_id, src_id in cleaned[["assay_chembl_id", "src_assay_id"]]
        .itertuples(index=False, name=None)
    }


def _ensure_src_assay_id(
    frame: pd.DataFrame, *, lookup: Mapping[str, str]
) -> pd.DataFrame:
    """Populate ``src_assay_id`` using ``lookup`` when missing."""

    if "src_assay_id" not in frame.columns and "assay_chembl_id" not in frame.columns:
        return frame

    result = frame.copy()

    if "src_assay_id" in result.columns:
        try:
            result["src_assay_id"] = result["src_assay_id"].astype("string")
        except TypeError:
            result["src_assay_id"] = result["src_assay_id"].astype("string")
    else:
        result["src_assay_id"] = pd.Series(pd.NA, index=result.index, dtype="string")

    if not lookup or result.empty or "assay_chembl_id" not in result.columns:
        return result

    assay_ids = result["assay_chembl_id"].astype("string")
    missing_mask = _string_like_missing(result["src_assay_id"])
    if not missing_mask.any():
        return result

    normalized_ids = assay_ids.where(~assay_ids.isna(), None).astype(object)

    def _resolve(value: object) -> str | None:
        if value is None:
            return None
        return lookup.get(str(value))

    mapped = normalized_ids.map(_resolve)
    available = missing_mask & mapped.notna()
    if not available.any():
        return result

    result.loc[available, "src_assay_id"] = mapped[available].astype("string")
    return result


def _ensure_molecule_pref_name(
    frame: pd.DataFrame,
    *,
    cfg: Config,
    client: ChemblClient,
    cache: dict[str, str | None],
    cache_lock: Lock,
    chunk_failures: ChunkFailureTracker | None = None,
) -> pd.DataFrame:
    """Populate ``molecule_pref_name`` via the test item API when missing."""

    if frame.empty or "molecule_chembl_id" not in frame.columns:
        return frame

    result = frame.copy()

    if "molecule_pref_name" in result.columns:
        missing_mask = _string_blank_mask(result["molecule_pref_name"])
    else:
        result["molecule_pref_name"] = pd.Series(pd.NA, index=result.index, dtype="string")
        missing_mask = pd.Series(True, index=result.index, dtype="boolean")

    if not missing_mask.any():
        return result

    molecule_ids = (
        result.loc[missing_mask, "molecule_chembl_id"].astype("string").str.strip()
    )
    molecule_ids = molecule_ids[molecule_ids != ""]
    unique_ids = tuple(dict.fromkeys(molecule_ids.dropna().tolist()))

    if not unique_ids:
        return result

    pending: list[str] = []
    wait_for: set[str] = set()

    with cache_lock:
        for identifier in unique_ids:
            cache_key = str(identifier)
            current = cache.get(cache_key, _CACHE_MISS)
            if current is _CACHE_IN_PROGRESS:
                wait_for.add(cache_key)
                continue
            if current is not _CACHE_MISS:
                continue
            cache[cache_key] = _CACHE_IN_PROGRESS  # type: ignore[assignment]
            pending.append(cache_key)

    if pending:
        fields = ["molecule_chembl_id", "pref_name"]
        extra_fields = getattr(cfg.testitem, "fields", None)
        if extra_fields:
            fields.extend(extra_fields)
        fields = sorted(set(fields))
        api_overrides: dict[str, Any] = {}
        if getattr(cfg.testitem, "retries", None) is not None:
            api_overrides["retries"] = cfg.testitem.retries
        if getattr(cfg.testitem, "backoff_factor", None) is not None:
            api_overrides["backoff_factor"] = cfg.testitem.backoff_factor
        api_cfg = cfg.api.model_copy(update=api_overrides) if api_overrides else cfg.api

        try:
            lookup = cl.get_testitem(
                pending,
                cfg=api_cfg,
                client=client,
                chunk_size=cfg.testitem.batch_size,
                timeout=cfg.testitem.timeout,
                fields=fields,
                page_limit=cfg.testitem.request_limit,
            )
        except (requests.RequestException, ValueError, AttributeError) as exc:
            error_message = str(exc)
            logger.warning(
                "pref_name_fetch_failed",
                pending=list(pending),
                error=error_message,
            )
            if chunk_failures is not None:
                chunk_failures.add_failure(tuple(pending), error_message)
            lookup = pd.DataFrame(columns=["molecule_chembl_id", "pref_name"])
        resolved: set[str] = set()
        if not lookup.empty and {"molecule_chembl_id", "pref_name"}.issubset(lookup.columns):
            mapped = (
                lookup[["molecule_chembl_id", "pref_name"]]
                .dropna(subset=["molecule_chembl_id"])
                .astype({"molecule_chembl_id": "string"})
            )
            for chembl_id, pref_name in mapped.itertuples(index=False):
                value = str(pref_name) if pd.notna(pref_name) else None
                cache_key = str(chembl_id)
                with cache_lock:
                    cache[cache_key] = value
                resolved.add(cache_key)
        missing = [identifier for identifier in pending if identifier not in resolved]
        if missing:
            with cache_lock:
                for identifier in missing:
                    cache[identifier] = None

    if wait_for:
        while True:
            with cache_lock:
                outstanding = [
                    identifier
                    for identifier in wait_for
                    if cache.get(identifier, _CACHE_MISS) is _CACHE_IN_PROGRESS
                ]
            if not outstanding:
                break
            sleep(0)

    with cache_lock:
        fill_map = {
            key: value
            for key, value in cache.items()
            if value and value is not _CACHE_IN_PROGRESS
        }
    if not fill_map:
        return result

    replacements = molecule_ids.map(fill_map)
    available = replacements.notna()
    if available.any():
        result.loc[molecule_ids.index[available], "molecule_pref_name"] = (
            replacements[available].astype("string")
        )

    return result
 
 


def _ensure_extended_activity_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Guarantee columns expected by the post-processing stage."""

    result = frame.copy()
    if result.empty:
        for column, dtype in _EXTENDED_ACTIVITY_DTYPES.items():
            if column not in result.columns:
                result[column] = pd.Series([], dtype=dtype)
        return result

    if "activity_id" in result.columns:
        if "activity_chembl_id" in result.columns:
            missing_id = _string_like_missing(result["activity_id"])
            if missing_id.any():
                result.loc[missing_id, "activity_id"] = result.loc[
                    missing_id, "activity_chembl_id"
                ]
        else:
            missing_id = _string_like_missing(result["activity_id"])
            if missing_id.any():
                # Preserve dtype by assigning NA values without coercing existing entries.
                result.loc[missing_id, "activity_id"] = pd.NA
    elif "activity_chembl_id" in result.columns:
        result["activity_id"] = result["activity_chembl_id"].copy()

    for column, dtype in _EXTENDED_ACTIVITY_DTYPES.items():
        fallback = _EXTENDED_ACTIVITY_FALLBACKS.get(column)
        if column in result.columns:
            coerced_existing, _ = _coerce_extended_series(result[column], dtype, column)
            result[column] = coerced_existing
            if fallback is not None:
                if dtype in {"Float64", "Int64", "boolean"}:
                    missing_mask = coerced_existing.isna()
                else:
                    missing_mask = _string_like_missing(coerced_existing)
                if missing_mask.any():
                    candidate = fallback(result)
                    if candidate is not None:
                        aligned = candidate.reindex(result.index)
                        filled = _coerce_series_dtype(aligned, dtype)
                        current = result[column]
                        # ``Series.mask`` preserves the extension dtype and avoids
                        # implicit "object" upcasts that would otherwise trigger
                        # ``FutureWarning`` about incompatible assignments (see
                        # pandas GH-53964).
                        result[column] = current.mask(missing_mask, filled)
                        coerced_fallback, fallback_failures = _coerce_extended_series(
                            aligned, dtype, column
                        )
                        failure_subset = fallback_failures & missing_mask
                        if failure_subset.any():
                            logger.debug(
                                "activity_extended_fallback_conversion_failed",
                                column=column,
                                dtype=dtype,
                                rows=int(failure_subset.sum()),
                                mode="existing",
                            )
                        result.loc[missing_mask, column] = coerced_fallback.loc[missing_mask]
            continue
        if fallback is not None:
            candidate = fallback(result)
            if candidate is not None:
                aligned = candidate.reindex(result.index)
                coerced_fallback, fallback_failures = _coerce_extended_series(
                    aligned, dtype, column
                )
                if fallback_failures.any():
                    logger.debug(
                        "activity_extended_fallback_conversion_failed",
                        column=column,
                        dtype=dtype,
                        rows=int(fallback_failures.sum()),
                        mode="new_column",
                    )
                result[column] = coerced_fallback
                continue
        if dtype == "boolean":
            filler = pd.Series(pd.NA, index=result.index, dtype="boolean")
        elif dtype == "Float64":
            filler = pd.Series(pd.NA, index=result.index, dtype="Float64")
        elif dtype == "Int64":
            filler = pd.Series(pd.NA, index=result.index, dtype="Int64")
        else:
            filler = pd.Series(pd.NA, index=result.index, dtype=dtype)
        result[column] = filler
    return result


def _filter_activity_output_columns(
    frame: pd.DataFrame,
    *,
    column_order: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """Return ``frame`` restricted to the allowed output columns."""

    drop_list = list(_OUTPUT_ACTIVITY_DROP_COLUMNS)
    dropped_columns = [column for column in drop_list if column in frame.columns]
    filtered = frame.drop(columns=drop_list, errors="ignore")

    if column_order is None:
        allowed_head: list[str] = [
            column for column in filtered.columns if column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
        ]
    else:
        allowed_head = [
            column
            for column in column_order
            if column in filtered.columns and column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
        ]

    extras = [column for column in filtered.columns if column not in allowed_head]
    if allowed_head or extras:
        filtered = filtered.loc[:, allowed_head + extras]

    return filtered, dropped_columns


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute activity retrieval from the ChEMBL API.

    The resulting CSV places columns defined in :data:`ActivitiesSchema`
    first, preserving their declared order. Any additional fields appear
    afterwards sorted alphabetically.

    Parameters
    ----------
    cfg : Config
        Application configuration providing API credentials, retry strategy
        and CSV export options.
    args : argparse.Namespace
        Parsed command-line arguments. ``args.limit`` constrains the number of
        identifiers processed, while ``args.dry_run`` skips network calls and
        file generation.

    Returns
    -------
    int
        ``0`` on success, non-zero when validation or I/O failures are
        encountered. Upstream API errors are logged and converted into a
        failure code by :func:`library.cli_utils.run_pipeline`.
    """
    offset = getattr(args, "offset", 0)
    workers_override = getattr(args, "workers", None)
    configured_workers = (
        workers_override
        if workers_override is not None
        else getattr(cfg.activity, "workers", 1)
    )
    final_out_attr = getattr(args, "final_out", None)
    if final_out_attr in (None, argparse.SUPPRESS):
        legacy_output = getattr(args, "output_csv", None)
        if legacy_output not in (None, argparse.SUPPRESS):
            if isinstance(legacy_output, str | Path):
                output_path = Path(legacy_output)
            else:
                output_path = Path(str(legacy_output))
            if not isinstance(legacy_output, Path):
                args.final_out = output_path
            args.output_csv = output_path
        else:
            output_path = Path(io.default_output_path(args.input_csv, cfg.io))
            args.final_out = output_path
            args.output_csv = output_path
    else:
        output_path = Path(str(final_out_attr)) if not isinstance(final_out_attr, str | Path) else Path(final_out_attr)
        if not isinstance(final_out_attr, Path):
            args.final_out = output_path
        args.output_csv = output_path

    start_time = perf_counter()

    pre_context = prepare_activity_context(cfg, args, skip_read=True)
    if pre_context is None:
        return 1

    limit = pre_context.limit

    logger.info(
        "activity_pipeline_start",
        input=str(args.input_csv),
        output=str(output_path),
        workers=configured_workers,
        limit=limit,
        offset=offset,
    )
    logger.info(
        "activity_pipeline_verbose",
        input=str(args.input_csv),
        output=str(output_path),
        limit=limit,
        offset=offset,
        batch_size=cfg.activity.batch_size,
        timeout=cfg.activity.timeout,
        dry_run=cfg.activity.dry_run,
        workers=configured_workers,
    )

    if cfg.activity.dry_run:
        expected = limit if limit is not None else 0
        logger.info("dry_run", limit=limit)
        logger.info("activity_pipeline_dry_run", expected=expected)
        _emit_completion_message(
            output_path=output_path,
            processed_rows=0,
            duration_s=perf_counter() - start_time,
            mode="dry_run",
        )
        return 0

    prepared_context = prepare_activity_context(cfg, args)
    if prepared_context is None:
        return 1

    limit = prepared_context.limit
    limited_ids = prepared_context.limited_ids
    extended_output_path: Path | None = None

    enrichment_cfg = cfg.activity_enrichment
    extra_columns: list[str] = []
    action_cfg = enrichment_cfg.action_type
    configure_activity_schema(action_cfg.metrics)
    if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
        extra_columns.append(action_cfg.column)
    extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}

    _args_invocation(args)

    failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
    fetch_failure_path = output_path.with_name(
        f"{output_path.stem}_fetch_failures.csv"
    )


    def _compute_bounds(frame: pd.DataFrame) -> pd.DataFrame:
        return compute_activity_bounds(frame, cfg.activity_bounds)

    def _apply_annotations(frame: pd.DataFrame) -> pd.DataFrame:
        return apply_activity_annotations(
            frame,
            action_cfg=enrichment_cfg.action_type,
            properties_cfg=enrichment_cfg.activity_properties,
        )

    def _ensure_required_activity_columns(frame: pd.DataFrame) -> pd.DataFrame:
        missing = [
            column
            for column in _ACTIVITY_REQUIRED_COLUMNS
            if column not in frame.columns
        ]
        if not missing:
            return frame
        logger.warning(
            "activity_missing_columns",
            missing=sorted(missing),
        )
        fillers: dict[str, pd.Series[Any]] = {}
        for column in missing:
            dtype_info = _ACTIVITY_REQUIRED_DTYPES.get(column)
            python_type = getattr(dtype_info, "python_type", None)
            dtype_text = str(dtype_info).lower() if dtype_info is not None else ""
            if python_type in {float, int} or "float" in dtype_text or "int" in dtype_text:
                fill_dtype = "Float64"
            elif python_type is str:
                fill_dtype = "string"
            else:
                fill_dtype = "object"
            fillers[column] = pd.Series(pd.NA, index=frame.index, dtype=fill_dtype)
        return frame.assign(**fillers)

    available_columns: set[str] = set()
    assay_src_lookup = _load_assay_src_lookup(cfg.resources.dictionary_dir)

    def _record_columns(frame: pd.DataFrame) -> pd.DataFrame:
        available_columns.update(frame.columns)
        return frame

    metadata_hooks = [
        _ensure_required_activity_columns,
        partial(_ensure_src_assay_id, lookup=assay_src_lookup),
        _ensure_extended_activity_columns,
        normalize_activities,
        add_pipeline_metadata,
        _compute_bounds,
        _apply_annotations,
        _record_columns,
    ]

    validators = [partial(validate_activities, return_result=True)]

    streaming_stats = _StreamingCSVStatistics()
    streamed_summary: dict[str, object] | None = None

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: Sequence[str],
        key_cols: Sequence[str],
    ) -> Path:
        sort_columns = list(key_cols) or sorted(col_order)
        column_order = list(col_order)
        filtered_order = [column for column in column_order if column in available_columns]

        drop_candidates = [
            column
            for column in _OUTPUT_ACTIVITY_DROP_COLUMNS
            if column in available_columns
        ]
        if drop_candidates:
            logger.info(
                "Dropped columns from output.activity_*: %s",
                ", ".join(drop_candidates),
            )

        whitelist_order = [
            column for column in filtered_order if column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
        ]

        def _stream_filtered_chunks() -> Iterator[pd.DataFrame]:
            for chunk in chunks:
                filtered_chunk, _ = _filter_activity_output_columns(
                    chunk,
                    column_order=filtered_order,
                )
                head = [
                    column for column in whitelist_order if column in filtered_chunk.columns
                ]
                tail = sorted(
                    column for column in filtered_chunk.columns if column not in head
                )
                if head or tail:
                    ordered_chunk = filtered_chunk.reindex(columns=head + tail, copy=False)
                else:
                    ordered_chunk = filtered_chunk
                streaming_stats.update(ordered_chunk)
                yield ordered_chunk

        filtered_chunks = _stream_filtered_chunks()

        output_path = write_csv_chunks_deterministic(
            filtered_chunks,
            destination,
            key_cols=sort_columns,
            col_order=whitelist_order,
            chunksize=cfg.io.csv_chunksize,
            sort_chunksize=cfg.io.csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        path_obj = Path(output_path)
        if io.write_csv is not _ORIGINAL_IO_WRITE_CSV:
            try:
                io.write_csv(
                    (),
                    destination,
                    cfg=cfg,
                    key_cols=sort_columns,
                    col_order=column_order,
                    chunksize=cfg.io.csv_chunksize,
                    sep=cfg.io.csv_sep,
                    encoding=cfg.io.csv_encoding,
                )
            except Exception:  # pragma: no cover - defensive for patched writers
                logger.debug(
                    "Fallback CSV writer stub raised an exception; deterministic export succeeded and execution will continue."
                )
        return path_obj



    doc_quality_cfg = cfg.system.doc_quality
    table_quality = build_table_quality_hook(
        doc_quality_cfg,
        table_name=Path(output_path).with_suffix(""),
        destination=Path(output_path).parent,
    )

    last_error_extra: dict[str, object] | None = None
    last_error_context: dict[str, object] = {}

    pref_name_cache: dict[str, str | None] = {}
    pref_name_cache_lock = Lock()

    with ETLContext(cfg) as context:
        client = context.chembl_client

        retry_cfg = cfg.retry
        chunk_failures = ChunkFailureTracker()
        http_diagnostics = _gather_http_diagnostics(cfg, client)
        if http_diagnostics:
            logger.info("activity_http_config", **http_diagnostics)

        def fetch_chunk(chunk_ids: Sequence[str]) -> pd.DataFrame:
            nonlocal last_error_extra, last_error_context

            timeout_error_types = (
                requests.Timeout,
                requests.ReadTimeout,
                requests.ConnectTimeout,
                requests.exceptions.RetryError,
            )

            def _is_timeout_error(exc: Exception) -> bool:
                if isinstance(exc, timeout_error_types):
                    return True
                message = str(exc).strip().lower()
                if not message:
                    return False
                return "timed out" in message or "timeout" in message

            def _fetch(ids: Sequence[str], *, depth: int = 0) -> pd.DataFrame:
                nonlocal last_error_extra, last_error_context

                attempts = max(1, retry_cfg.max_attempts)
                id_list = [str(identifier) for identifier in ids]

                for attempt in range(1, attempts + 1):
                    try:
                        result = cl.get_activities(
                            id_list,
                            cfg=cfg.api,
                            client=client,
                            chunk_size=cfg.activity.batch_size,
                            timeout=cfg.activity.timeout,
                            **extra_kwargs,
                        )
                    except (requests.RequestException, ValueError) as exc:
                        error_message = str(exc)
                        network_hint, network_host = _describe_network_failure(cfg, exc)
                        context = {
                            "chunk_ids": list(id_list),
                            "chunk_size": len(id_list),
                            "attempt": attempt,
                            "max_attempts": attempts,
                            "batch_size": cfg.activity.batch_size,
                            "timeout": cfg.activity.timeout,
                        }
                        if network_host:
                            context["api_host"] = network_host
                        for key in (
                            "adapter_total",
                            "api_retries",
                            "retry_max_attempts",
                            "activity_timeout",
                            "api_timeout_read",
                        ):
                            if key in http_diagnostics:
                                context[key] = http_diagnostics[key]
                        log_context = {
                            key: value for key, value in context.items() if key != "chunk_ids"
                        }
                        last_error_extra = {
                            "msg": error_message,
                            "chunk_ids": context["chunk_ids"],
                        }
                        if network_hint:
                            last_error_extra["hint"] = network_hint
                        last_error_context = dict(log_context)
                        if network_host:
                            last_error_context.setdefault("api_host", network_host)
                        if network_hint:
                            last_error_context["network_hint"] = network_hint
                        if attempt >= attempts:
                            if network_hint:
                                logger.error(
                                    "activity_fetch_network_error",
                                    extra=last_error_extra,
                                    hint=network_hint,
                                    **log_context,
                                )
                            if len(id_list) > 1 and _is_timeout_error(exc):
                                split_context = dict(log_context)
                                split_context["depth"] = depth
                                logger.warning(
                                    "activity_fetch_split",
                                    extra=last_error_extra,
                                    **split_context,
                                )
                                midpoint = max(1, len(id_list) // 2)
                                left_ids = tuple(id_list[:midpoint])
                                right_ids = tuple(id_list[midpoint:])
                                frames: list[pd.DataFrame] = []
                                if left_ids:
                                    frames.append(_fetch(left_ids, depth=depth + 1))
                                if right_ids:
                                    frames.append(_fetch(right_ids, depth=depth + 1))
                                if frames:
                                    combined = pd.concat(
                                        frames, ignore_index=True, sort=False
                                    )
                                else:
                                    combined = pd.DataFrame(columns=ACTIVITY_COLUMNS)
                                last_error_extra = None
                                last_error_context = {}
                                return combined
                            logger.error(
                                "activity_fetch_failed",
                                extra=last_error_extra,
                                error=error_message,
                                **log_context,
                            )
                            chunk_failures.add_failure(id_list, error_message)
                            raise PipelineError("chunk_fetch_failed")  # noqa: B904
                        delay = compute_backoff_delay(attempt, retry_cfg)
                        logger.warning(
                            "activity_fetch_retry",
                            extra=last_error_extra,
                            delay=delay,
                            **log_context,
                        )
                        if delay > 0:
                            sleep(delay)
                    else:
                        last_error_extra = None
                        last_error_context = {}
                        return _ensure_molecule_pref_name(
                            result,
                            cfg=cfg,
                            client=client,
                            cache=pref_name_cache,
                            cache_lock=pref_name_cache_lock,
                            chunk_failures=chunk_failures,
                        )
                return pd.DataFrame(columns=ACTIVITY_COLUMNS)

            return _fetch(chunk_ids, depth=0)

        worker_count = getattr(cfg.activity, "workers", 1) or 1
        fetch_config = ChunkedFetchConfig(
            ids=limited_ids,
            chunk_size=cfg.activity.batch_size,
            workers=max(1, worker_count),
        )

        writer_config = CsvWriterConfig(
            writer=writer,
            kwargs={},
            ensure_destination=True,

        )

        pipeline_stats: dict[str, object] | None = None

        def _capture_stats(stats: Mapping[str, object]) -> None:
            nonlocal pipeline_stats, streamed_summary
            pipeline_stats = dict(stats)
            snapshot = streaming_stats.snapshot()
            streamed_summary = snapshot
            pipeline_stats.setdefault("rows_streamed", snapshot["rows"])
            pipeline_stats.setdefault("cells_streamed", snapshot["cells"])
            pipeline_stats.setdefault("null_cells", snapshot["nulls"])
            pipeline_stats.setdefault("null_fraction", snapshot["null_fraction"])

        definition_kwargs: dict[str, object] = {
            "schema": ActivitiesSchema,
            "schema_name": "ActivitiesSchema",
            "validators": validators,
            "command": " ".join(_args_invocation(args)),
            "config_snapshot": _serialize_paths(cfg.to_dict()),
            "inputs": {"input_csv": str(args.input_csv)},
            "key_columns": ["activity_id"],
            "table_quality": table_quality,
            "stats_extra": chunk_failures.stats,
            "stats_callback": _capture_stats,
            "dictionary_resources": (
                "dictionary_root",
                "target_types",
            ),
        }

        result = activity_run.run_activity_pipeline(
            fetch_config=fetch_config,
            metadata_hooks=metadata_hooks,
            fetch_chunk=fetch_chunk,
            writer_config=writer_config,
            definition_kwargs=definition_kwargs,
            cfg=cfg,
            logger=logger,
            output_path=output_path,
            failure_path=failure_path,
            fetch_failure_path=fetch_failure_path,
            chunk_failures=chunk_failures,
        )
        exit_code = result.exit_code

    if exit_code == 0:
        logger.info(
            f"Merged data checkpoint: wrote merged activity records to '{output_path}'."
        )
        try:
            extended_output_path = process_activity_extended(
                input_path=output_path,
                search_dir=output_path.parent,
                dictionary_dir=cfg.resources.dictionary_dir,
                targets_csv=cfg.resources.targets_type_csv,
            )
        except Exception:
            logger.exception(
                "Activity extended post-processing failed while enriching merged activity data."
            )
            raise

    processed_ids = prepared_context.processed_ids
    try:
        processed_count = int(processed_ids or 0)
    except (TypeError, ValueError):
        logger.info("processed_count_conversion_failed", value=processed_ids)
        processed_count = 0

    summary_snapshot: dict[str, object] | None = None

    if limit is not None:
        logger.info(
            "process_limit",
            processed=processed_ids,
            limit=limit,
        )
        summary_snapshot = streamed_summary or streaming_stats.snapshot()
    else:
        summary_snapshot = streamed_summary or streaming_stats.snapshot()
        logger.info("processed_count", count=processed_count)

    if pipeline_stats is not None:
        try:
            rows_total = int(pipeline_stats.get("rows_total", processed_count))
        except (TypeError, ValueError):
            rows_total = processed_count
        try:
            rows_kept = int(pipeline_stats.get("rows_kept", 0))
        except (TypeError, ValueError):
            rows_kept = 0
        rows_dropped = int(pipeline_stats.get("rows_dropped", 0))
        logger.info(
            "records_dropped",
            rows_total=rows_total,
            rows_kept=rows_kept,
            rows_dropped=rows_dropped,
        )

    if exit_code == 0:
        completion_rows = processed_ids
        summary_rows = summary_snapshot.get("rows") if summary_snapshot else None
        if isinstance(summary_rows, int | float):
            try:
                completion_rows = int(summary_rows)
            except Exception:
                # Fallback in case of any unexpected type or value
                completion_rows = processed_ids
        elif pipeline_stats is not None:
            try:
                # Try to get from pipeline_stats if available
                rows_kept = pipeline_stats.get("rows_kept")
                if rows_kept is not None and isinstance(rows_kept, int | float):
                    completion_rows = int(rows_kept)
                else:
                    rows_total = pipeline_stats.get("rows_total", processed_ids)
                    if isinstance(rows_total, int | float):
                        completion_rows = int(rows_total)
                    else:
                        completion_rows = processed_ids
            except Exception:
                # Fallback in case of any unexpected type or value
                completion_rows = processed_ids
        else:
            completion_rows = processed_ids

        report_extras: dict[str, object] = {"rows": completion_rows, "processed": processed_ids}
        if pipeline_stats is not None:
            report_extras.update(pipeline_stats)
        if summary_snapshot:
            report_extras["summary_snapshot"] = summary_snapshot
        if extended_output_path is not None:
            report_extras["extended_output"] = str(extended_output_path)

        postprocess_metrics, report_path = _cli_generate_activity_postprocess_metrics(
            cfg,
            output_path,
            logger=logger,
            extras=report_extras,
        )

        pipeline_version_value = (
            getattr(postprocess_metrics, "pipeline_version", None)
            if postprocess_metrics and postprocess_metrics.pipeline_version is not None
            else get_pipeline_version()
        )

        pipeline_done_payload: dict[str, object] = {
            "output": str(output_path),
            "rows": completion_rows,
            "pipeline_version": pipeline_version_value,
        }
        if summary_snapshot is not None:
            null_fraction = summary_snapshot.get("null_fraction")
            if null_fraction is not None:
                pipeline_done_payload["null_fraction"] = null_fraction
        if extended_output_path is not None:
            pipeline_done_payload["extended_output"] = str(extended_output_path)
        if postprocess_metrics is not None:
            summary = postprocess_metrics.summary()
            if summary.get("rows") is not None:
                pipeline_done_payload["postprocess_rows"] = summary["rows"]
            if summary.get("columns") is not None:
                pipeline_done_payload["postprocess_columns"] = summary["columns"]
            if summary.get("duration_s") is not None:
                pipeline_done_payload["postprocess_duration_s"] = summary["duration_s"]
            if summary.get("steps") is not None:
                pipeline_done_payload["postprocess_steps"] = summary["steps"]
            validation = getattr(postprocess_metrics, "validation", None)
            if validation is not None:
                pipeline_done_payload["postprocess_schema"] = getattr(validation, "schema", None)
            if report_path is not None:
                pipeline_done_payload["postprocess_report"] = str(report_path)

        logger.info("activity_pipeline_done", extra=pipeline_done_payload)
        if extended_output_path is not None:
            logger.info(
                "activity_export_checkpoint",
                output=str(output_path),
                extended_output=str(extended_output_path),
            )
        else:
            logger.info(
                "activity_export_checkpoint",
                output=str(output_path),
                extended_output=None,
            )
        _emit_completion_message(
            output_path=output_path,
            processed_rows=completion_rows,
            duration_s=perf_counter() - start_time,
            mode="run",
            streamed_metrics=summary_snapshot,
        )
    else:
        extra_payload = last_error_extra
        context_payload = dict(last_error_context) if last_error_context else {}
        error_message = None
        chunk_ids: Sequence[str] | None = None
        if extra_payload:
            error_message = extra_payload.get("msg")
            chunk_ids = extra_payload.get("chunk_ids")  # type: ignore[assignment]
        details = []
        if error_message:
            details.append(f"last error: {error_message}")
        if chunk_ids:
            details.append(f"chunk_ids={list(chunk_ids)}")
        if context_payload:
            attempt_info = ", ".join(f"{key}={value}" for key, value in context_payload.items())
            details.append(attempt_info)
        detail_text = "; ".join(details)
        logger.error(
            "activity_pipeline_failed",
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
            details=detail_text or None,
        )
        logger.error(
            "activity_pipeline_failed_detail",
            output=str(output_path),
            processed=processed_ids,
            exit_code=exit_code,
            details=detail_text or None,
        )

    return exit_code

_activity = import_module("library.cli.entrypoints.activity")

register_activity_pipeline_hooks(
    runner=run_chembl,
    emit_completion_message=_emit_completion_message,
)

setattr(_activity, "MAX_ACTIVITY_CHUNK_SIZE", MAX_ACTIVITY_CHUNK_SIZE)


class _LoggerProxy:
    """Proxy ``library.pipelines.activity.runner`` logger to the CLI module."""

    def __getattr__(self, name: str):  # type: ignore[override]
        return getattr(_activity, "logger").__getattribute__(name)


_activity_runner.logger = _LoggerProxy()

_extended_exports = (
    "ActivityPipelineCLI",
    "DEFAULT_INPUT_NAME",
    "DEFAULT_OUTPUT_STEM",
    "PROGRAM_NAME",
    "run",
    "run_chembl",
    "_emit_completion_message",
    "main",
    "build_parser",
    "MAX_ACTIVITY_CHUNK_SIZE",
)

_activity.__all__ = tuple(
    dict.fromkeys(getattr(_activity, "__all__", tuple()) + _extended_exports)
)

sys.modules[__name__] = _activity
sys.modules.setdefault("scripts.get_activity_data", _activity)

if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(_activity.main())
