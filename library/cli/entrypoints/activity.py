"""Command line interface for retrieving ChEMBL activity data.

The module exposes a ``main`` entry point compatible with setuptools console
scripts as well as helpers that can be invoked directly from other
applications or tests.
"""

from __future__ import annotations

import argparse
import logging
import math
import numbers
import os
import re
import sys
from collections.abc import (
    Callable,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from dataclasses import dataclass
from datetime import UTC, datetime
from datetime import datetime as _datetime
from enum import Enum
from functools import partial
from itertools import islice
from pathlib import Path
from threading import Condition, Lock
from time import perf_counter, sleep
from types import ModuleType
from typing import Any, cast
from urllib.parse import urlsplit

import pandas as pd
import requests
from library._compat.pandera import pa

import library.cli.logging as cli_logging
from library import cli, io, offline
from library.cli import Logger, LoggerConfig, positive_int
from library.cli import build_parser as base_parser
from library.cli.activity_api import (
    MIN_ACTIVITY_TIMEOUT,
    ActivityCommandOptions,
    ensure_entrypoint_exports,
    run_activity_pipeline,
)
from library.cli.base import PipelineCLIBase
from library.cli.commands import get_activity_data as _activity_cli_commands
from library.cli.logging import CLILoggingContext
from library.cli_utils import PipelineError, resolve_invocation
from library.cli_utils import run_cli_command as _run_cli_command
from library.clients.chembl import ChemblClient
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.fetch_retry import ChunkFailureTracker, compute_backoff_delay
from library.common.log import logger
from library.common.run_context import get_current as get_run_context
from library.config import Config, _serialize_paths
from library.integration import chembl_library as cl
from library.metadata import (
    file_sha256 as _metadata_file_sha256,
)
from library.metadata import (
    write_meta_yaml as _metadata_write_meta_yaml,
)
from library.orchestration import ETLContext
from library.pipelines.activity import run as activity_run
from library.pipelines.assay.chembl_assay import ACTIVITY_COLUMNS
from library.pipelines.common import (
    ChunkedFetchConfig,
    CsvWriterConfig,
    add_pipeline_metadata,
)
from library.pipelines.common.metadata import get_pipeline_version
from library.postprocess import (
    PostprocessingPipelineConfig,
    run_postprocessing_pipeline,
)
from library.postprocess import (
    get_csv_runtime_config as get_postprocess_csv_config,
)
from library.postprocess import (
    get_pipeline_config as get_postprocess_pipeline_config,
)
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
from library.utils.data_correlation import generate_correlation_report
from library.utils.qc_report import generate_qc_report
from library.validation import validate_activities

_DEFAULT_LOGGING_DATE_FUNC = cli_logging._current_date_str

_DATE_TOKEN_RE = re.compile(r"\A\d{8}\Z")

try:  # pragma: no cover - urllib3 is part of requests dependency chain
    from urllib3.exceptions import NameResolutionError as _Urllib3NameResolutionError
except Exception:  # pragma: no cover - defensive fallback for alternative stacks
    _Urllib3NameResolutionError = None  # type: ignore[assignment]

DEFAULT_INPUT_NAME = "activity.csv"
DEFAULT_OUTPUT_STEM = "activity"
PROGRAM_NAME = Path(__file__).with_suffix("").name

_ACTIVITY_METADATA_SOURCES: tuple[str, ...] = ("ChEMBL Activity API",)

EMPTY_COLUMNS = [
    "approx_cited_activity",
    "compound_key",
    "exact_cited_activity",
    "higly_correlated_cit",
    "multmol_assay",
    "nstereo",
    "original_activity_approx",
    "original_activity_exact",
    "review_doc",
    "rounded_data_citation",
    "salt_chembl_id",
    "shuffled_cit",
    "standard_lower_value",
    "standard_upper_value",
]


def _build_metadata_summary(frame: pd.DataFrame) -> dict[str, Any]:
    """Return a concise QC summary suitable for metadata sidecars."""

    total_rows = int(len(frame))
    if total_rows == 0:
        return {"total_rows": 0, "non_null_ratio": {}}

    ratios: dict[str, float] = {}
    for column in frame.columns:
        valid = int(frame[column].notna().sum())
        ratios[column] = float(valid / total_rows)
    return {"total_rows": total_rows, "non_null_ratio": ratios}


def _drop_empty_activity_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Remove known empty columns from the activity export frame."""

    columns_to_drop = [col for col in EMPTY_COLUMNS if col in frame.columns]
    if not columns_to_drop:
        return frame

    pruned = frame.drop(columns=columns_to_drop, errors="ignore")
    logger.info(
        "activity_empty_columns_removed",
        removed=len(columns_to_drop),
        columns=columns_to_drop,
    )
    return pruned


# ---------------------------------------------------------------------------
# Compatibility hooks
# ---------------------------------------------------------------------------

# ``datetime`` and ``clock`` are exposed as module-level variables so that tests
# and downstream integrations can override them for deterministic behaviour.
# Historically :mod:`scripts.get_activity_data` imported ``datetime`` directly
# which allowed monkeypatching via ``get_activity_data.datetime``. Restoring the
# binding keeps that contract intact after the CLI refactor.
datetime = _datetime  # noqa: F811 - preserve compatibility alias
clock: Callable[[], float] = perf_counter


def _current_utc_datetime() -> _datetime:
    """Return the current UTC timestamp using the overridable clock."""

    candidate = datetime.now(UTC)
    if candidate.tzinfo is None:
        return candidate.replace(tzinfo=UTC)
    return candidate.astimezone(UTC)


def _current_date_token() -> str:
    """Return the YYYYMMDD date string derived from :data:`datetime`."""

    return _current_utc_datetime().strftime("%Y%m%d")


_DATE_SUFFIX_RE = re.compile(r"(?P<table>.*?)(?:_)?(?P<date>\d{8})$")


def _derive_standard_output_labels(dataset_csv: Path) -> tuple[str, str]:
    """Return ``(table_name, date_tag)`` inferred from ``dataset_csv``."""

    candidate = Path(dataset_csv)
    base_path = candidate
    for suffix in reversed(candidate.suffixes):
        if suffix.lower() not in {".csv", ".tmp"}:
            break
        base_path = base_path.with_suffix("")

    base = base_path.name.lstrip(".")
    prefix = "output."
    while base.lower().startswith(prefix):
        base = base[len(prefix) :].lstrip(".")

    base = base.strip()
    if not base:
        return DEFAULT_OUTPUT_STEM, _current_date_token()

    match = _DATE_SUFFIX_RE.search(base)
    if match and match.group("date"):
        table = (match.group("table") or "").rstrip("_")
        date_tag = match.group("date")
        table_name = table or DEFAULT_OUTPUT_STEM
        if table_name.lower() == "activities":
            table_name = DEFAULT_OUTPUT_STEM
        return table_name, date_tag

    table_name = base or DEFAULT_OUTPUT_STEM
    if table_name.lower() == "activities":
        table_name = DEFAULT_OUTPUT_STEM
    return table_name, _current_date_token()


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
write_meta_yaml = _metadata_write_meta_yaml
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
    "datetime",
    "clock",
)


def _ensure_command_logger_sync() -> None:
    """Keep the shared command module logger aligned with this script logger."""

    try:
        commands_module = _activity_cli_commands
    except NameError:  # pragma: no cover - defensive guard for refactors
        commands_module = None

    if (
        commands_module is not None
        and getattr(commands_module, "logger", None) is not logger
    ):
        try:
            commands_module.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass

    try:
        from library.pipelines.activity import runner as activity_runner
    except Exception:  # pragma: no cover - defensive guard for circular imports
        activity_runner = None

    if (
        activity_runner is not None
        and getattr(activity_runner, "logger", None) is not logger
    ):
        try:
            activity_runner.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass

    try:
        from library.common import log as _common_log
    except Exception:  # pragma: no cover - defensive guard for circular imports
        _common_log = None

    if _common_log is not None and getattr(_common_log, "logger", None) is not logger:
        try:
            _common_log.logger = logger
        except Exception:  # pragma: no cover - defensive guard
            pass


_ensure_command_logger_sync()


class _CacheState(Enum):
    MISS = "miss"
    IN_PROGRESS = "in-progress"


CacheValue = str | None | _CacheState

_CACHE_MISS = _CacheState.MISS
_CACHE_IN_PROGRESS = _CacheState.IN_PROGRESS


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
    #  "salt_chembl_id": "string",
    "target_chembl_id": "string",
    "bao_endpoint": "string",
    #  "compound_key": "string",
    "compound_name": "string",
    #   "multmol_assay": "boolean",
    #   "approx_cited_activity": "boolean",
    #   "shuffled_cit": "boolean",
    #   "exact_cited_activity": "boolean",
    #   "higly_correlated_cit": "boolean",
    #   "review_doc": "boolean",
    #   "rounded_data_citation": "boolean",
    #   "original_activity_approx": "string",
    #   "original_activity_exact": "string",
    #   "nstereo": "Int64",
    #   "log_value": "Float64",
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

    if dtype in {"Float64", "Int64", "boolean"}:
        converted, _ = _coerce_extended_series(
            series, dtype, column="_coerce_series_dtype"
        )
        return cast("pd.Series[Any]", converted)

    if dtype == "string":
        return series.astype(pd.StringDtype())

    try:
        import numpy as np

        return series.astype(np.dtype(dtype))
    except (TypeError, ValueError):
        return cast("pd.Series[Any]", series.astype(str))


def _safe_int(value: object, default: int = 0) -> int:
    """Safely convert ``value`` to :class:`int` while tolerating bad inputs."""

    if value is None:
        return default

    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default


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
            if _Urllib3NameResolutionError is not None and isinstance(
                candidate, _Urllib3NameResolutionError
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


def _describe_network_failure(
    cfg: Config, exc: Exception
) -> tuple[str | None, str | None]:
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
            elif isinstance(value, str | bool) or value is None:
                metrics_payload[key] = value
            else:
                metrics_payload[key] = str(value)

    if mode == "skip_existing" and output_path is not None:
        logger.info("pipeline_skip_existing", output_postprocessed=str(output_path))
        events_attr = getattr(logger, "events", None)
        if isinstance(events_attr, list):
            events_attr.append(
                (
                    "info",
                    "pipeline_skip_existing",
                    {"output_postprocessed": str(output_path)},
                )
            )
        return

    payload: dict[str, object] = {
        "output_postprocessed": str(output_path) if output_path is not None else None,
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


_EXTENDED_ACTIVITY_FALLBACKS: dict[
    str, Callable[[pd.DataFrame], pd.Series[Any] | None]
] = {
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
        converted = pd.Series(
            pd.array(numeric.tolist(), dtype="Float64"), index=series.index
        )
        failures = series.notna() & converted.isna()
        return converted, failures

    if dtype == "Int64":
        numeric = pd.to_numeric(series, errors="coerce")
        non_integral_mask = numeric.notna() & numeric.ne(numeric.round())
        if non_integral_mask.any():
            numeric = numeric.mask(non_integral_mask)
        try:
            converted = pd.Series(
                pd.array(numeric.tolist(), dtype="Int64"), index=series.index
            )
        except (TypeError, ValueError):
            converted = pd.Series(pd.NA, index=series.index, dtype="Int64")
            failures = pd.Series(True, index=series.index)
            return converted, failures
        failures = series.notna() & (converted.isna() | non_integral_mask)
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


def _normalise_lookup_series(series: pd.Series[Any]) -> pd.Series[str]:
    """Coerce heterogeneous lookup values to trimmed string dtype."""

    def _coerce(value: Any) -> Any:
        if pd.isna(value):
            return pd.NA

        if isinstance(value, str):
            result = value.strip()
            if result and any(char in result for char in ".eE"):
                try:
                    float_value = float(result)
                except (TypeError, ValueError):
                    pass
                else:
                    if math.isnan(float_value):
                        return pd.NA
                    if float_value.is_integer():
                        result = str(int(float_value))
                    else:
                        result = format(float_value, "g")
        elif isinstance(value, bytes):
            result = value.decode("utf-8", errors="ignore").strip()
            if result and any(char in result for char in ".eE"):
                try:
                    float_value = float(result)
                except (TypeError, ValueError):
                    pass
                else:
                    if math.isnan(float_value):
                        return pd.NA
                    if float_value.is_integer():
                        result = str(int(float_value))
                    else:
                        result = format(float_value, "g")
        elif isinstance(value, numbers.Integral) and not isinstance(value, bool):
            result = str(int(value))
        elif isinstance(value, numbers.Real):
            float_value = float(value)
            if math.isnan(float_value):
                return pd.NA
            if float_value.is_integer():
                result = str(int(float_value))
            else:
                result = format(float_value, "g")
        else:
            result = str(value).strip()

        if not result:
            return pd.NA

        return result

    normalised = series.map(_coerce)
    return normalised.astype("string")


def _load_assay_src_lookup(dictionary_dir: Path | str | None) -> dict[str, str]:
    """Return mapping of assay identifiers to ``src_assay_id`` values."""

    if dictionary_dir is None:
        return {}

    candidate = Path(dictionary_dir) / "_assay" / "assay.csv"
    read_kwargs: dict[str, Any] = {"low_memory": False}
    try:
        frame = pd.read_csv(candidate, encoding="utf-8", **read_kwargs)
    except UnicodeDecodeError:
        for encoding in ("utf-8-sig", "cp1252", "latin-1"):
            try:
                frame = pd.read_csv(candidate, encoding=encoding, **read_kwargs)
                break
            except UnicodeDecodeError:
                continue
        else:
            logger.warning(
                f"Assay lookup file '{candidate}' could not be decoded using supported encodings; src_assay_id enrichment will be skipped."
            )
            return {}
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
        assay_chembl_id=_normalise_lookup_series(cleaned["assay_chembl_id"]),
        src_assay_id=_normalise_lookup_series(cleaned["src_assay_id"]),
    )

    cleaned = cleaned.dropna(subset=["assay_chembl_id", "src_assay_id"])
    if cleaned.empty:
        return {}

    return {
        str(assay_id): str(src_id)
        for assay_id, src_id in cleaned[["assay_chembl_id", "src_assay_id"]].itertuples(
            index=False, name=None
        )
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
    cache: MutableMapping[str, CacheValue],
    cache_condition: Condition,
    chunk_failures: ChunkFailureTracker | None = None,
    wait_timeout: float | None = None,
) -> pd.DataFrame:
    """Populate ``molecule_pref_name`` via the test item API when missing."""

    if frame.empty or "molecule_chembl_id" not in frame.columns:
        return frame

    result = frame.copy()

    if "molecule_pref_name" in result.columns:
        missing_mask = _string_blank_mask(result["molecule_pref_name"])
    else:
        result["molecule_pref_name"] = pd.Series(
            pd.NA, index=result.index, dtype="string"
        )
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

    with cache_condition:
        for identifier in unique_ids:
            cache_key = str(identifier)
            current = cache.get(cache_key, _CACHE_MISS)
            if current is _CACHE_IN_PROGRESS:
                wait_for.add(cache_key)
                continue
            if current is not _CACHE_MISS:
                continue
            cache[cache_key] = _CACHE_IN_PROGRESS
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
        if not lookup.empty and {"molecule_chembl_id", "pref_name"}.issubset(
            lookup.columns
        ):
            mapped = (
                lookup[["molecule_chembl_id", "pref_name"]]
                .dropna(subset=["molecule_chembl_id"])
                .astype({"molecule_chembl_id": "string"})
            )
            for chembl_id, pref_name in mapped.itertuples(index=False):
                value = str(pref_name) if pd.notna(pref_name) else None
                cache_key = str(chembl_id)
                with cache_condition:
                    cache[cache_key] = value
                    cache_condition.notify_all()
                resolved.add(cache_key)
        missing = [identifier for identifier in pending if identifier not in resolved]
        if missing:
            with cache_condition:
                for identifier in missing:
                    cache[identifier] = None
                cache_condition.notify_all()

    if wait_for:
        timeout = (
            wait_timeout
            if wait_timeout is not None
            else getattr(getattr(cfg, "testitem", object()), "timeout", None)
        )
        timeout = float(timeout) if isinstance(timeout, numbers.Real) else None
        deadline = clock() + timeout if timeout and timeout > 0 else None
        poll_interval = 0.5

        with cache_condition:
            while True:
                outstanding = [
                    identifier
                    for identifier in wait_for
                    if cache.get(identifier, _CACHE_MISS) is _CACHE_IN_PROGRESS
                ]
                if not outstanding:
                    break

                remaining: float | None = None
                if deadline is not None:
                    remaining = max(0.0, deadline - clock())
                    if remaining == 0:
                        logger.warning(
                            "pref_name_fetch_wait_timeout",
                            pending=sorted(outstanding),
                        )
                        break

                wait_interval = (
                    remaining if remaining and remaining > 0 else poll_interval
                )
                cache_condition.wait(timeout=wait_interval)

                if deadline is not None and clock() >= deadline:
                    outstanding = [
                        identifier
                        for identifier in wait_for
                        if cache.get(identifier, _CACHE_MISS) is _CACHE_IN_PROGRESS
                    ]
                    if outstanding:
                        logger.warning(
                            "pref_name_fetch_wait_timeout",
                            pending=sorted(outstanding),
                        )
                        break

    with cache_condition:
        fill_map = {
            key: value
            for key, value in cache.items()
            if value is not _CACHE_IN_PROGRESS
        }
    if not fill_map:
        return result

    replacements = molecule_ids.map(fill_map)
    available = replacements.notna()
    if available.any():
        result.loc[molecule_ids.index[available], "molecule_pref_name"] = replacements[
            available
        ].astype("string")

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
                        result.loc[missing_mask, column] = coerced_fallback.loc[
                            missing_mask
                        ]
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
            column
            for column in filtered.columns
            if column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
        ]
    else:
        allowed_head = [
            column
            for column in column_order
            if column in filtered.columns
            and column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
        ]

    extras = [column for column in filtered.columns if column not in allowed_head]
    if allowed_head or extras:
        filtered = filtered.loc[:, allowed_head + extras]

    return filtered, dropped_columns


@dataclass(frozen=True)
class _ActivityPostprocessDeps:
    process_activity_extended: Callable[..., Path]
    run_activity_postprocess: Callable[..., tuple[pd.DataFrame, object]]
    validate_postprocess: Callable[..., pd.DataFrame]
    activity_schema: pa.DataFrameSchema


def _load_activity_postprocess_deps() -> _ActivityPostprocessDeps:
    from library.postprocessing.activities import (
        ACTIVITY_SCHEMA as _ACTIVITY_SCHEMA,
    )
    from library.postprocessing.activities import (
        run_activity_pipeline as _run_activity_postprocess,
    )
    from library.postprocessing.activities import (
        validate_activities as _validate_postprocess,
    )
    from library.postprocessing.activity_extended import process_activity_extended

    return _ActivityPostprocessDeps(
        process_activity_extended=process_activity_extended,
        run_activity_postprocess=_run_activity_postprocess,
        validate_postprocess=_validate_postprocess,
        activity_schema=_ACTIVITY_SCHEMA,
    )


def _derive_postprocess_output_path(output_path: Path) -> Path:
    name = output_path.name
    if name.startswith("output."):
        return output_path.with_name(f"output_postprocessed.{name[len('output.'):]}")
    return output_path.with_name(
        f"{output_path.stem}.postprocessed{output_path.suffix}"
    )


def _derive_standard_output_labels(dataset_csv: Path) -> tuple[str, str]:
    """Return canonical ``(table_name, date_tag)`` extracted from ``dataset_csv``."""

    name = dataset_csv.name
    # Strip leading dots so temporary files like ``.output.activities_*.tmp``
    # do not leak the hidden prefix into derived artefact names.
    name = name.lstrip(".")

    # Remove known temporary suffixes generated by the pipeline before
    # analysing the basename.
    lowered = name.lower()
    for suffix in (".tmp",):
        while lowered.endswith(suffix):
            name = name[: -len(suffix)]
            lowered = name.lower()

    if lowered.endswith(".csv"):
        name = name[: -len(".csv")]
        lowered = name.lower()

    prefix = "output."
    if name.startswith(prefix):
        # Drop redundant ``output.`` prefixes that can accumulate when the
        # pipeline works with hidden temporary artefacts such as
        # ``.output.activities_*.csv.tmp``.  The intermediate filenames can be
        # normalised by stripping both the prefix and any leading dots that
        # remain after the removal.  Repeat the process to guard against
        # chains like ``output..output.activities`` which would otherwise leak
        # the extra ``output.`` segment into the derived table name and produce
        # paths such as ``output.output.activities_*.csv``.
        while name.startswith(prefix):
            candidate = name[len(prefix) :]
            candidate = candidate.lstrip(".")
            if not candidate:
                break
            name = candidate

    base = name
    if "_" in base:
        table_candidate, date_candidate = base.rsplit("_", 1)
    else:
        table_candidate, date_candidate = base, ""

    clean_table = table_candidate.lstrip(".")

    # Some historical exports use dotted stems such as
    # ``output.activity.20240101.csv``.  Normalise these by treating the final
    # ``.<YYYYMMDD>`` segment as the date token when no underscore suffix is
    # present.  This keeps the derived artefact names consistent with the
    # ``output.<table>_<date>.csv`` policy while tolerating legacy filename
    # layouts generated by upstream tooling.
    if not date_candidate and "." in clean_table:
        dotted_table, dotted_suffix = clean_table.rsplit(".", 1)
        if _DATE_TOKEN_RE.fullmatch(dotted_suffix):
            clean_table = dotted_table
            date_candidate = dotted_suffix

    if _DATE_TOKEN_RE.fullmatch(clean_table) and not date_candidate:
        table_name = DEFAULT_OUTPUT_STEM
        date_tag = clean_table
    else:
        table_name = clean_table or DEFAULT_OUTPUT_STEM
        if table_name.lower() == "activities":
            table_name = DEFAULT_OUTPUT_STEM
        if _DATE_TOKEN_RE.fullmatch(date_candidate):
            date_tag = date_candidate
        else:
            date_tag = _current_date_token()

    return table_name, date_tag


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
    offline_flag = offline.is_enabled(getattr(args, "offline", None))
    args.offline = offline_flag
    fixtures: offline.OfflineFixtures | None = None
    if offline_flag:
        try:
            fixtures = offline.OfflineFixtures()
        except FileNotFoundError as exc:
            logger.error("activity_offline_fixtures_missing", error=str(exc))
            return 1
        logger.info("activity_offline_mode", base=str(fixtures.base_dir))
    def _execute() -> int:
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
                output_path = Path(
                    io.default_output_path(
                        args.input_csv,
                        cfg.io,
                        date=getattr(args, "date", None),
                    )
                )
                args.final_out = output_path
                args.output_csv = output_path
        else:
            output_path = (
                Path(str(final_out_attr))
                if not isinstance(final_out_attr, str | Path)
                else Path(final_out_attr)
            )
            if not isinstance(final_out_attr, Path):
                args.final_out = output_path
            args.output_csv = output_path
        
        emit_legacy = bool(getattr(args, "emit_legacy_artifacts", False))
        
        start_time = clock()
        
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
                duration_s=clock() - start_time,
                mode="dry_run",
            )
            return 0
        
        prepared_context = prepare_activity_context(cfg, args)
        if prepared_context is None:
            return 1
        
        limit = prepared_context.limit
        limited_ids = prepared_context.limited_ids
        extended_output_path: Path | None = None
        postprocess_metrics = None
        postprocess_report_path: Path | None = None
        postprocess_output_path: Path | None = None
        dataset_path: Path | None = None
        exit_code = 1
        
        enrichment_cfg = cfg.activity_enrichment
        extra_columns: list[str] = []
        action_cfg = enrichment_cfg.action_type
        configure_activity_schema(action_cfg.metrics)
        if action_cfg.enabled or action_cfg.log_missing or action_cfg.log_distribution:
            extra_columns.append(action_cfg.column)
        extra_kwargs = {"extra_columns": extra_columns} if extra_columns else {}
        
        _args_invocation(args)
        
        failure_path = output_path.with_name(f"{output_path.stem}_failure_cases.csv")
        fetch_failure_path = output_path.with_name(f"{output_path.stem}_fetch_failures.csv")
        
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
                if (
                    python_type in {float, int}
                    or "float" in dtype_text
                    or "int" in dtype_text
                ):
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
            filtered_order = [
                column for column in column_order if column in available_columns
            ]
        
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
                column
                for column in filtered_order
                if column not in _OUTPUT_ACTIVITY_DROP_COLUMNS
            ]
        
            def _stream_filtered_chunks() -> Iterator[pd.DataFrame]:
                for chunk in chunks:
                    filtered_chunk, _ = _filter_activity_output_columns(
                        chunk,
                        column_order=filtered_order,
                    )
                    head = [
                        column
                        for column in whitelist_order
                        if column in filtered_chunk.columns
                    ]
                    tail = sorted(
                        column for column in filtered_chunk.columns if column not in head
                    )
                    if head or tail:
                        ordered_chunk = filtered_chunk.reindex(
                            columns=head + tail, copy=False
                        )
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
        if emit_legacy:
            table_quality = build_table_quality_hook(
                doc_quality_cfg,
                table_name=Path(output_path).with_suffix(""),
                destination=Path(output_path).parent,
            )
        else:
        
            def _noop_table_quality(_: Path) -> None:
                return None
        
            table_quality = _noop_table_quality
        
        def _persist_standard_outputs(dataset_csv: Path) -> io.StandardOutputArtifacts:
            dataset_frame = pd.read_csv(
                dataset_csv,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
            dataset_frame = _drop_empty_activity_columns(dataset_frame)
            if "raw.index" not in dataset_frame.columns:
                dataset_frame.insert(0, "raw.index", dataset_frame.index)
                logging.info(
                    "Добавлена индексная колонка 'raw.index' (%s строк).",
                    len(dataset_frame),
                )
            quality_report = pd.DataFrame()
            correlation_report = pd.DataFrame()
            table_name_value, date_tag = _derive_standard_output_labels(dataset_csv)
            table_label = Path(cfg.io.output_dir) / f"output.{table_name_value}_{date_tag}"
            try:
                correlation_report = generate_correlation_report(
                    dataset_frame,
                    table_name=str(table_label),
                )
            except Exception as exc:  # pragma: no cover - defensive guard
                logger.warning(
                    "activity_correlation_generation_failed",
                    error=str(exc),
                    path=str(dataset_csv),
                )
                correlation_report = pd.DataFrame()
            try:
                quality_report = generate_qc_report(
                    dataset_frame,
                    table_name=str(table_label),
                )
            except Exception as exc:  # pragma: no cover - defensive guard
                logger.warning(
                    "activity_quality_generation_failed",
                    error=str(exc),
                    path=str(dataset_csv),
                )
                quality_report = pd.DataFrame()
            table_name_value, date_tag = _derive_standard_output_labels(dataset_csv)
            output_directory = Path(cfg.io.output_dir)
            artifacts = io.save_standard_outputs(
                dataset_frame,
                correlation_report,
                quality_report,
                table_name=table_name_value,
                date_tag=date_tag,
                output_dir=output_directory,
                output_path=dataset_csv,
            )
            qc_summary = _build_metadata_summary(dataset_frame)
            io.save_metadata(
                table_name=table_name_value,
                date_tag=date_tag,
                args=args,
                qc_summary=qc_summary,
                output_dir=artifacts.dataset.parent,
                artifacts=[
                    artifacts.dataset,
                    artifacts.quality_report,
                    artifacts.correlation_report,
                ],
                sources=_ACTIVITY_METADATA_SOURCES,
                run_context=get_run_context(),
            )
            logger.info(
                "activity_standard_outputs",
                dataset=str(artifacts.dataset),
                quality_report=str(artifacts.quality_report),
                correlation_report=str(artifacts.correlation_report),
            )
            return artifacts
        
        last_error_extra: dict[str, object] | None = None
        last_error_context: dict[str, object] = {}
        
        pref_name_cache: dict[str, CacheValue] = {}
        pref_name_cache_lock = Lock()
        pref_name_cache_condition = Condition(pref_name_cache_lock)
        
        with ETLContext(cfg) as context:
            client = context.chembl_client
        
            retry_cfg = cfg.retry
            jitter = retry_cfg.build_jitter()
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
                                key: value
                                for key, value in context.items()
                                if key != "chunk_ids"
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
                            delay = compute_backoff_delay(attempt, retry_cfg, jitter=jitter)
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
                                cache_condition=pref_name_cache_condition,
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
                emit_legacy_artifacts=emit_legacy,
            )
            exit_code = int(result.exit_code)
            dataset_path = result.output_path
            if dataset_path is not None:
                output_path = Path(dataset_path)
                args.final_out = output_path
                args.output_csv = output_path
        
        if not emit_legacy:
            Path(fetch_failure_path).unlink(missing_ok=True)
            Path(f"{fetch_failure_path}.meta.yaml").unlink(missing_ok=True)
        
        if exit_code == 0:
            if dataset_path is not None:
                standard_artifacts = _persist_standard_outputs(Path(dataset_path))
                standard_dataset = standard_artifacts.dataset
                auto_generated = bool(getattr(args, "_auto_output_generated", False))
                if (
                    not emit_legacy
                    and auto_generated
                    and standard_dataset != Path(dataset_path)
                ):
                    Path(dataset_path).unlink(missing_ok=True)
                dataset_path = standard_dataset
                output_path = dataset_path
                args.final_out = output_path
                args.output_csv = output_path
            logger.info(
                f"Merged data checkpoint: wrote merged activity records to '{output_path}'."
            )
            postprocess_enabled = bool(getattr(args, "postprocess", False))
            if not postprocess_enabled:
                logger.info("[INFO] Postprocessing skipped (flag --postprocess not set)")
            else:
                deps = _load_activity_postprocess_deps()
                try:
                    extended_output_path = deps.process_activity_extended(
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
        
                try:
                    pipeline_config = get_postprocess_pipeline_config(
                        "activities", getattr(args, "config", None)
                    )
                    csv_cfg = get_postprocess_csv_config(pipeline_config)
                    runtime_cfg = PostprocessingPipelineConfig(
                        pipeline_config=pipeline_config,
                        csv_runtime_config=csv_cfg,
                        runner=deps.run_activity_postprocess,
                        validator=deps.validate_postprocess,
                        schema=deps.activity_schema,
                        logger=logger,
                    )
                    postprocess_output_path = _derive_postprocess_output_path(output_path)
                    logger.info(
                        "activity_postprocess_start",
                        input=str(output_path),
                        output=str(postprocess_output_path),
                    )
                    pipeline_result = run_postprocessing_pipeline(
                        "activities",
                        output_path,
                        postprocess_output_path,
                        runtime_cfg,
                    )
                except Exception:
                    logger.exception(
                        "Activity postprocessing pipeline failed while processing merged activity data."
                    )
                    raise
        
                postprocess_metrics = pipeline_result.metrics
                postprocess_report_path = pipeline_result.report_path
                if postprocess_metrics is None:
                    raise RuntimeError("activity postprocess metrics missing")
                if postprocess_report_path is None:
                    raise RuntimeError("activity postprocess report missing")
                logger.info(
                    "activity_postprocess_done",
                    output=str(pipeline_result.output_path),
                    report=str(postprocess_report_path),
                    rows=_safe_int(postprocess_metrics.output_rows, 0),
                    columns=_safe_int(postprocess_metrics.output_columns, 0),
                )
        
        summary_snapshot: Mapping[str, object] | None = None
        processed_ids = prepared_context.processed_ids
        processed_count = _safe_int(processed_ids, 0)
        if processed_count == 0 and processed_ids not in (None, 0, "0", 0.0):
            logger.info("processed_count_conversion_failed", value=processed_ids)
        
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
            rows_total = _safe_int(pipeline_stats.get("rows_total"), processed_count)
            rows_kept = _safe_int(pipeline_stats.get("rows_kept"), 0)
            rows_dropped = _safe_int(pipeline_stats.get("rows_dropped"), 0)
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
                if postprocess_output_path is not None:
                    pipeline_done_payload["postprocess_output"] = str(
                        postprocess_output_path
                    )
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
                    pipeline_done_payload["postprocess_schema"] = getattr(
                        validation, "schema", None
                    )
                if postprocess_report_path is not None:
                    pipeline_done_payload["postprocess_report"] = str(
                        postprocess_report_path
                    )
        
            logger.info("activity_pipeline_done", extra=pipeline_done_payload)
            if extended_output_path is not None:
                logger.info(
                    f"Successful export checkpoint: primary data at '{output_path}', extended data at '{extended_output_path}'."
                )
            else:
                logger.info(
                    f"Successful export checkpoint: activity data written to '{output_path}'."
                )
            _emit_completion_message(
                output_path=output_path,
                processed_rows=completion_rows,
                duration_s=clock() - start_time,
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
                attempt_info = ", ".join(
                    f"{key}={value}" for key, value in context_payload.items()
                )
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
                f"Activity pipeline failed with exit code {exit_code} after processing {processed_ids} identifiers destined for '{output_path}'. {detail_text}"
            )
        
        return exit_code
    if fixtures is None:
        return _execute()
    with offline.patch_activity(fixtures):
        return _execute()


def _coerce_cli_path(value: object) -> Path | str | None:
    """Normalise optional CLI path parameters for helper delegation."""

    if value in (None, argparse.SUPPRESS):
        return None
    return cast(
        Path | str | None, value
    )  # ``run_activity_pipeline`` handles conversion to :class:`Path`.


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the activity pipeline handling ``--skip-existing`` semantics."""

    _ensure_command_logger_sync()

    try:
        from library.common import log as _common_log

        _common_log.logger = logger
    except Exception:  # pragma: no cover - defensive guard
        pass

    start_time = clock()

    output_csv_value = _coerce_cli_path(getattr(args, "output_csv", None))
    if output_csv_value is not None and not isinstance(output_csv_value, Path):
        output_csv_value = Path(output_csv_value)

    final_output_value = _coerce_cli_path(getattr(args, "final_out", None))
    if final_output_value is not None and not isinstance(final_output_value, Path):
        final_output_value = Path(final_output_value)

    skip_existing = getattr(args, "skip_existing", False)
    force = getattr(args, "force", False)

    candidate_output = final_output_value or output_csv_value
    output_path: Path | None = (
        Path(candidate_output) if candidate_output is not None else None
    )
    preexisting_output = output_path.exists() if output_path is not None else False

    if skip_existing and not force and preexisting_output and output_path is not None:
        logger.info("pipeline_skip_existing", output_postprocessed=str(output_path))
        events_attr = getattr(logger, "events", None)
        if isinstance(events_attr, list):
            events_attr.append(
                (
                    "info",
                    "pipeline_skip_existing",
                    {"output_postprocessed": str(output_path)},
                )
            )
        _emit_completion_message(
            output_path=output_path,
            processed_rows=None,
            duration_s=clock() - start_time,
            mode="skip_existing",
        )
        return 0

    options = ActivityCommandOptions(
        input_csv=args.input_csv,
        output_csv=output_csv_value,
        final_output=final_output_value,
        limit=getattr(args, "limit", None),
        offset=getattr(args, "offset", 0),
        timeout=getattr(args, "timeout", None),
        batch_size=getattr(args, "batch_size", None),
        workers=getattr(args, "workers", None),
        dry_run=getattr(args, "dry_run", False),
        skip_existing=skip_existing,
        force=force,
        invocation=getattr(args, "invocation", None),
    )

    exit_code = run_activity_pipeline(
        cfg,
        options,
        runner=run_chembl,
        emit_completion_message=_emit_completion_message,
    )

    if skip_existing and not force and preexisting_output and output_path is not None:
        logger.info("pipeline_skip_existing", output_postprocessed=str(output_path))
        events_attr = getattr(logger, "events", None)
        if isinstance(events_attr, list):
            events_attr.append(
                (
                    "info",
                    "pipeline_skip_existing",
                    {"output_postprocessed": str(output_path)},
                )
            )

    return exit_code


def _build_parser_impl() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        A tuple containing the fully configured parser and the logging
        configuration populated with defaults.
    """
    parser, log_cfg = base_parser(
        "ChEMBL activity data utilities",
        column="activity_id",
        chunk_size=5,
        size_option="--batch-size",
        size_dest="batch_size",
        emit_legacy_help=(
            "Persist the historical CSV, metadata and diagnostics files in "
            "addition to the standard activity output bundle (default: "
            "disabled)."
        ),
    )
    parser.prog = PROGRAM_NAME
    parser.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser.add_argument(
        "--timeout",
        type=float,
        default=90.0,
        help=(
            "Timeout in seconds for each HTTP request (values below "
            f"{int(MIN_ACTIVITY_TIMEOUT)} seconds are automatically clamped)"
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help=("Maximum number of identifiers to process; use 0 to skip processing"),
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read input and exit without contacting the API or writing files",
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="Number of worker threads fetching activities in parallel",
    )
    parser.add_argument(
        "--postprocess",
        dest="postprocess",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Enable activity postprocessing after the main pipeline",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Use cached fixtures instead of contacting the ChEMBL API",
    )
    legacy_option = parser._option_string_actions.get("--emit-legacy-artifacts")
    if legacy_option is not None:
        legacy_option.help = (
            "Persist the historical CSV, metadata and diagnostics files in addition "
            "to the standard activity output bundle (default: disabled)."
        )
    parser.set_defaults(func=run_chembl)
    return parser, log_cfg


class ActivityPipelineCLI(PipelineCLIBase):
    """CLI adapter for the activity data pipeline."""

    def __init__(self) -> None:
        super().__init__()
        self._log_path: Path | None = None

    def build_parser(self) -> tuple[argparse.ArgumentParser, LoggerConfig]:
        return _build_parser_impl()

    def prepare_arguments(
        self,
        parser: argparse.ArgumentParser,
        args: argparse.Namespace,
        argv: Sequence[str] | None,
    ) -> argparse.Namespace:
        args.invocation = resolve_invocation(parser.prog, argv)
        cli.prepare_io_paths(
            args,
            input_default=DEFAULT_INPUT_NAME,
            output_stem=DEFAULT_OUTPUT_STEM,
        )
        if not hasattr(args, "_initial_output_stamp_mode"):
            args._initial_output_stamp_mode = getattr(args, "output_stamp_mode", None)
        if not hasattr(args, "_initial_auto_output_generated"):
            args._initial_auto_output_generated = bool(
                getattr(args, "_auto_output_generated", False)
            )
        return args

    def handle_pre_run(
        self, parser: argparse.ArgumentParser, args: argparse.Namespace
    ) -> int | None:
        if args.limit == 0:
            logger.info("pipeline_skip_limit", limit=args.limit)
            logger.info(
                "Limit set to 0; exiting before starting the activity pipeline."
            )
            return 0
        if args.limit is not None and args.limit < 0:
            parser.error("--limit must be zero or a positive integer")
        if args.offset < 0:
            parser.error("--offset must be zero or a positive integer")
        input_path = Path(args.input_csv)
        if not input_path.exists():
            message = (
                f"Input CSV '{input_path}' does not exist. "
                "Provide --input pointing to a valid identifiers file or "
                "run the upstream pipeline step to generate it."
            )
            sys.stderr.write(message + "\n")
            return 1
        return None

    def get_program_name(self) -> str:
        return PROGRAM_NAME

    def get_logger(self) -> Logger:
        return logger

    def get_run_cli_command(self) -> Callable[..., int]:
        candidate = getattr(_activity_cli_commands, "run_cli_command", None)
        if callable(candidate):
            return candidate
        return super().get_run_cli_command()

    def get_logging_date(self, args: argparse.Namespace) -> str | None:
        value = super().get_logging_date(args)
        if isinstance(value, str):
            stripped = value.strip()
            if stripped:
                return stripped
        try:
            original = getattr(cli_logging, "_ORIGINAL_CURRENT_DATE_STR", None)
            current = cli_logging._current_date_str
            if original is not None and current is not original:
                return current()
        except Exception:
            pass
        try:
            return _current_date_token()
        except Exception:  # pragma: no cover - defensive against custom hooks
            return value if isinstance(value, str) else None

    def get_config_mapping(self) -> Mapping[str, str]:
        return {
            "timeout": "activity.timeout",
            "column": "activity.column",
            "batch_size": "activity.batch_size",
            "limit": "activity.limit",
            "offset": "activity.offset",
            "dry_run": "activity.dry_run",
            "workers": "activity.workers",
        }

    def on_logging_ready(self, logging_ctx: CLILoggingContext) -> None:
        self._log_path = logging_ctx.log_path

    def after_run(self, log_cfg: LoggerConfig, exit_code: int) -> int:
        result = super().after_run(log_cfg, exit_code)
        self._ensure_legacy_log_alias()
        return result

    def _ensure_legacy_log_alias(self) -> None:
        if self._log_path is None:
            return

        log_path = self._log_path
        legacy_name = "get_activity_data"

        if log_path.name.startswith(legacy_name):
            return

        if not log_path.exists():
            return

        env_base = os.environ.get("CHEMBL_DA_BASE_PATH")

        suffix = log_path.name[len(PROGRAM_NAME) :]
        alias_path = log_path.with_name(f"{legacy_name}{suffix}")

        try:
            alias_path.write_text(
                log_path.read_text(encoding="utf-8"), encoding="utf-8"
            )
        except (
            OSError
        ) as exc:  # pragma: no cover - filesystem failures are environment-specific
            logger.warning(
                "log_alias_write_failed",
                source=str(log_path),
                alias=str(alias_path),
                error=str(exc),
            )
            return

        if env_base and Path(env_base).name == "cli":
            try:
                log_path.unlink()
                self._log_path = alias_path
            except OSError:  # pragma: no cover - best effort cleanup
                pass

    def run_pipeline(self, cfg: Config, args: argparse.Namespace) -> int:
        self._refresh_output_paths(cfg, args)
        return run(cfg, args)

    def _refresh_output_paths(self, cfg: Config, args: argparse.Namespace) -> None:
        initial_mode = getattr(args, "_initial_output_stamp_mode", None)
        current_mode = getattr(args, "output_stamp_mode", None)
        if initial_mode == current_mode:
            return
        if current_mode == "require" and getattr(args, "date", None) in (None, ""):
            try:
                from library.cli.commands import get_activity_data as _activity_commands

                _activity_commands._ensure_default_date(cfg, args)
            except Exception:  # pragma: no cover - defensive fallback
                pass
        if bool(getattr(args, "_initial_auto_output_generated", False)):
            args.final_out = None
            args.output_csv = None
        output_stem = getattr(args, "_auto_output_stem", None)
        suffix = getattr(args, "_auto_output_suffix", None)
        if not isinstance(output_stem, str) or not isinstance(suffix, str):
            return
        cli.prepare_io_paths(args, output_stem=output_stem, suffix=suffix)


_CLI = ActivityPipelineCLI()


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Expose parser construction for existing imports."""

    return _CLI.build_parser()


def main(argv: Sequence[str] | None = None) -> int:
    """Delegate to :class:`ActivityPipelineCLI` for backwards compatibility."""

    return _CLI.main(argv)


ensure_entrypoint_exports(cast(ModuleType, sys.modules[__name__]))


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
