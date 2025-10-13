"""Reusable components for the test item data acquisition pipeline."""

from __future__ import annotations

import os
import sys
import tempfile
import threading
import time
import traceback
import weakref
from collections import OrderedDict, deque
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import chain, islice, tee
from pathlib import Path
from types import MappingProxyType
from typing import IO, Any

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import io
from library.clients import pubchem as pc
from library.common.log import logger
from library.common.metadata import (
    Stats,
    file_sha256,
    record_quality_failure,
    write_meta_yaml,
)
from library.common.run_context import get_current
from library.common.sidecar import SidecarErrors, resolve_failure_chunk_size
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    PubChemCfg,
    RetryCfg,
    TestitemBatchRetryCfg,
    TestitemMoleculeEnrichmentCfg,
    _serialize_paths,
)
from library.integration.chembl_client import ChemblClient
from library.orchestration import ETLContext
from library.pipelines.assay.chembl_assay import TESTITEM_STRUCTURE_COLUMNS
from library.qa.reporting import build_table_quality_hook
from library.qa.validation import validate_testitems
from library.schemas import TestitemsSchema, normalize_testitems
from library.table_quality import _apply_sampling_and_filters
from library.utils import data_correlation, qc_report

from .catalog import (
    PARENT_LOOKUP_SOURCE_SKIPPED,
    ParentLookupStats,
    _merge_parent_stats,
    prepare_parent_enrichment,
    run_parent_enrichment,
)

try:  # pragma: no cover - exercised during package bootstrap
    from . import testitem_enrichment
except ImportError:  # pragma: no cover - fallback for partial initialisation
    from . import enrichment as testitem_enrichment


@lru_cache(maxsize=1)
def _chembl_library() -> Any:
    """Return the lazily imported ChEMBL integration module."""

    from library.integration import chembl_library

    return chembl_library


_FETCH_ERROR_SAMPLE_SIZE = 10
# ``_MISSING_IDENTIFIER_STATS_SAMPLE_SIZE`` controls how many missing identifiers
# are persisted in execution stats to avoid generating extremely large log
# payloads or metadata entries for massive gaps in the source dataset.
_MISSING_IDENTIFIER_LOG_SAMPLE_SIZE = 10
_MISSING_IDENTIFIER_STATS_SAMPLE_SIZE = 100
# ``finalize_output`` may consume placeholder rows from streaming sources, so
# keep the chunk size moderate to avoid large peak allocations when many IDs are
# missing from the source dataset.
MISSING_IDENTIFIER_PLACEHOLDER_CHUNK_SIZE = 1000
_DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE = 10
_PLACEHOLDER_CONTACT_EMAIL = "contact@example.org"
_DEFAULT_TABLE_NAME = "testitem"
RAW_INDEX_COLUMN = "raw.index"

_REQUESTED_IDENTIFIER_LOG_SAMPLE_SIZE = 10
_REQUESTED_IDENTIFIER_DUPLICATE_SAMPLE_SIZE = 10
_REQUESTED_IDENTIFIER_DUPLICATE_TRACK_LIMIT = 2048


class RequestedIdsSnapshot(Sequence[str]):
    """Persist requested identifiers without holding them all in memory."""

    def __init__(
        self,
        *,
        log_sample_size: int = _REQUESTED_IDENTIFIER_LOG_SAMPLE_SIZE,
        duplicate_sample_size: int = _REQUESTED_IDENTIFIER_DUPLICATE_SAMPLE_SIZE,
        duplicate_track_limit: int = _REQUESTED_IDENTIFIER_DUPLICATE_TRACK_LIMIT,
    ) -> None:
        self._log_sample_size = max(0, int(log_sample_size))
        self._duplicate_sample_size = max(0, int(duplicate_sample_size))
        self._duplicate_track_limit = max(0, int(duplicate_track_limit))
        self._log_sample: deque[str] = deque()
        self._duplicate_sample: list[str] = []
        self._duplicate_count = 0
        self._seen: OrderedDict[str, None] = OrderedDict()
        self._count = 0
        self._writer: IO[str] | None = None
        self._path: Path | None = None
        self._finalizer: weakref.finalize | None = None
        self._finished = False
        self._unique_original: OrderedDict[str, str] = OrderedDict()

    def __len__(self) -> int:
        return self._count

    def __iter__(self) -> Iterator[str]:
        return self._iterate_all()

    def __getitem__(self, index: int | slice) -> str | list[str]:
        if isinstance(index, slice):
            start, stop, step = index.indices(self._count)
            if step == 1 and start == 0 and stop == self._count:
                return list(self._iterate_all())
            return [self[i] for i in range(start, stop, step)]
        if index < 0:
            index += self._count
        if index < 0 or index >= self._count:
            raise IndexError("requested identifier index out of range")
        for idx, value in enumerate(self._iterate_all()):
            if idx == index:
                return value
        raise IndexError("requested identifier index out of range")

    def append(self, identifier: str | None) -> None:
        text_identifier = str(identifier)
        if self._log_sample_size and len(self._log_sample) < self._log_sample_size:
            self._log_sample.append(text_identifier)

        normalized = text_identifier.strip().upper()
        if normalized:
            if normalized in self._seen:
                self._duplicate_count += 1
                if (
                    self._duplicate_sample_size
                    and len(self._duplicate_sample) < self._duplicate_sample_size
                ):
                    self._duplicate_sample.append(text_identifier)
            else:
                self._seen[normalized] = None
                self._seen.move_to_end(normalized)
                if (
                    self._duplicate_track_limit
                    and len(self._seen) > self._duplicate_track_limit
                ):
                    self._seen.popitem(last=False)
                if normalized not in self._unique_original:
                    self._unique_original[normalized] = text_identifier

        self._ensure_writer()
        if self._writer is None:
            return
        self._writer.write(text_identifier + "\n")
        self._count += 1

    def finish(self) -> None:
        if self._finished:
            return
        self._finished = True
        if self._writer is not None:
            self._writer.flush()
            self._writer.close()

    def iter_for_fetch(self) -> Iterator[str]:
        self.finish()
        return self._read_from_disk()

    @property
    def sample(self) -> tuple[str, ...]:
        return tuple(self._log_sample)

    @property
    def duplicate_count(self) -> int:
        return self._duplicate_count

    @property
    def duplicate_sample(self) -> tuple[str, ...]:
        return tuple(self._duplicate_sample)

    @property
    def tracked_unique_count(self) -> int:
        return len(self._seen)

    @property
    def unique_map(self) -> Mapping[str, str]:
        return MappingProxyType(self._unique_original)

    def _ensure_writer(self) -> None:
        if self._writer is not None or self._finished:
            return
        temp = tempfile.NamedTemporaryFile(
            mode="w+",
            encoding="utf-8",
            newline="\n",
            delete=False,
            prefix="chembl_requested_ids_",
        )
        self._writer = temp
        self._path = Path(temp.name)
        self._finalizer = weakref.finalize(self, self._cleanup_path, self._path)

    def _read_from_disk(self) -> Iterator[str]:
        if self._path is None:
            return iter(())

        def _reader() -> Iterator[str]:
            if self._path is None:
                return
            with self._path.open("r", encoding="utf-8") as handle:
                for line in handle:
                    yield line.rstrip("\n")

        return _reader()

    def _iterate_all(self) -> Iterator[str]:
        self.finish()
        return self._read_from_disk()

    @staticmethod
    def _cleanup_path(path: Path) -> None:
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def ensure_raw_index_column(frame: pd.DataFrame) -> bool:
    """Ensure ``frame`` contains the raw index column as the leading int64 field."""

    if RAW_INDEX_COLUMN in frame.columns:
        frame.loc[:, RAW_INDEX_COLUMN] = frame[RAW_INDEX_COLUMN].astype(
            "int64", copy=False
        )
        return False

    raw_index = pd.Series(
        frame.index.to_numpy(dtype="int64", copy=False),
        index=frame.index,
        name=RAW_INDEX_COLUMN,
        dtype="int64",
    )
    frame.insert(0, RAW_INDEX_COLUMN, raw_index)
    return True


_DEFAULT_METADATA_SOURCES = (
    "ChEMBL Molecule API",
    "PubChem Compound API",
)

_QC_RATIO_COLUMNS = (
    "pubchem_cid",
    "pubchem_canonical_smiles",
)


@dataclass(frozen=True)
class PubChemAugmentationContext:
    """Capture the parameters required for a fallback PubChem enrichment run."""

    pubchem_cfg: PubChemCfg
    api_cfg: ApiCfg
    retry_cfg: RetryCfg
    client: ChemblClient
    timeout: float
    fields: Sequence[str] | None
    request_limit: int


def _normalise_output_labels(
    output: Path | str,
    *,
    default_table: str = _DEFAULT_TABLE_NAME,
    fallback_date: str | None = None,
) -> tuple[str, str]:
    """Return canonical table name and date tag derived from ``output``.

    Intermediate filenames produced by the pipeline (for example
    ``.output.testitem_20240101.csv.tmp``) are normalised so that downstream
    artefacts re-use the canonical ``testitem`` label. When ``fallback_date`` is
    omitted the helper mirrors the previous behaviour by falling back to the
    current UTC date.
    """

    table_name, date_tag = io.derive_output_labels(
        output,
        default_table=default_table,
        fallback_date=fallback_date,
    )
    if table_name.lower() == "testitems":
        table_name = default_table
    return table_name, date_tag


# Canonical dataset stem used for standard outputs.
_DEFAULT_OUTPUT_TABLE = "testitem"

_TABLE_NAME_ALIASES: dict[str, str] = {"testitems": _DEFAULT_OUTPUT_TABLE}

_PUBCHEM_OPTIONAL_COLUMNS = frozenset(
    {
        "pubchem_canonical_smiles",
        "pubchem_cid",
        "pubchem_inchi",
        "pubchem_inchikey",
        "pubchem_isomeric_smiles",
        "pubchem_iupac_name",
        "pubchem_molecular_formula",
    }
)
_PUBCHEM_KEY_COLUMNS = frozenset({"pubchem_cid"})
_SALT_OPTIONAL_COLUMN = "salt_chembl_id"
_DEFAULT_TABLE_NAME = "testitem"


@dataclass
class ReadInputIdsResult:
    """Container holding the identifier iterator and a diagnostic sample."""

    ids_iter: Iterator[str]
    sample_ids: tuple[str, ...]


@dataclass(frozen=True)
class TestitemPipelineOptions:
    """Parameters controlling CLI execution of the test item pipeline."""

    input_csv: Path
    output_csv: Path | None = None
    limit: int | None = None
    offset: int | None = None
    emit_legacy_artifacts: bool = False
    pubchem_enabled: bool | None = None
    date: str | None = None


def _resolve_metadata_sources(pubchem_enabled: bool | None) -> list[str]:
    """Return data provenance sources for the metadata payload."""

    sources = [str(_DEFAULT_METADATA_SOURCES[0])]
    if pubchem_enabled:
        sources.append(str(_DEFAULT_METADATA_SOURCES[1]))
    return sources


def _extract_metadata_parameters(
    cfg: Config,
    options: TestitemPipelineOptions | None,
    *,
    emit_legacy_artifacts: bool,
    pubchem_enabled: bool | None,
) -> dict[str, object]:
    """Collect pipeline parameters for metadata emission."""

    parameters: dict[str, object] = {
        "emit_legacy_artifacts": bool(emit_legacy_artifacts),
    }

    def _coerce_optional(value: object | None) -> object | None:
        return value if value is not None else None

    limit = _coerce_optional(getattr(options, "limit", None))
    if limit is None:
        limit = _coerce_optional(getattr(getattr(cfg, "testitem", None), "limit", None))
    if limit is not None:
        parameters["limit"] = int(limit)

    offset = _coerce_optional(getattr(options, "offset", None))
    if offset is None:
        offset = _coerce_optional(getattr(getattr(cfg, "testitem", None), "offset", None))
    if offset is not None:
        parameters["offset"] = int(offset)

    retry_cfg = getattr(cfg, "retry", None)
    retries = getattr(retry_cfg, "max_attempts", None) if retry_cfg is not None else None
    if retries is None:
        batch_retry = getattr(getattr(cfg, "testitem", None), "batch_retry", None)
        if batch_retry is not None:
            retries = getattr(batch_retry, "max_attempts", None)
    if retries is not None:
        parameters["retries"] = int(retries)

    pubchem_toggle = getattr(options, "pubchem_enabled", None)
    if pubchem_toggle is not None:
        parameters["pubchem_enabled"] = bool(pubchem_toggle)
    else:
        cfg_pubchem = getattr(cfg, "pubchem", None)
        enabled = getattr(cfg_pubchem, "enable", None) if cfg_pubchem is not None else None
        if enabled is not None:
            parameters["pubchem_enabled"] = bool(enabled)
        elif pubchem_enabled is not None:
            parameters["pubchem_enabled"] = bool(pubchem_enabled)

    return parameters


def _build_qc_summary(frame: pd.DataFrame) -> dict[str, Any]:
    """Return a compact QC summary for the metadata payload."""

    total_rows = int(frame.shape[0])
    summary: dict[str, Any] = {"total_rows": total_rows}
    if total_rows <= 0:
        return summary

    ratios: dict[str, float] = {}
    for column in _QC_RATIO_COLUMNS:
        if column not in frame.columns:
            continue
        non_null = frame[column].notna().sum()
        ratios[column] = round(float(non_null) / float(total_rows), 4)
    if ratios:
        summary["non_null_ratio"] = ratios
    return summary


def _metadata_generated_at(candidate: str | None = None) -> str | None:
    """Return the preferred ``generated_at`` timestamp for metadata files."""

    if candidate:
        stripped = str(candidate).strip()
        if stripped:
            return stripped

    context = get_current()
    if context is not None:
        generated_at = getattr(context, "generated_at", None)
        if isinstance(generated_at, str):
            stripped = generated_at.strip()
            if stripped:
                return stripped
    return None


def _write_primary_metadata(
    *,
    dataset_frame: pd.DataFrame | None,
    dataset_path: Path,
    quality_path: Path,
    correlation_path: Path,
    table_name: str,
    date_tag: str,
    parameters: Mapping[str, object] | None,
    stats: Mapping[str, object] | None,
    pubchem_enabled: bool | None,
    qc_summary: Mapping[str, object] | None = None,
    generated_at: str | None = None,
) -> Path:
    """Persist the standard metadata sidecar next to ``dataset_path``."""

    outputs = [dataset_path.name, quality_path.name, correlation_path.name]
    canonical = table_name.lower()
    resolved_table = _TABLE_NAME_ALIASES.get(canonical, table_name)
    if resolved_table.lower() != _DEFAULT_OUTPUT_TABLE:
        resolved_table = _DEFAULT_OUTPUT_TABLE

    qc_payload = (
        _build_qc_summary(dataset_frame)
        if dataset_frame is not None
        else dict(qc_summary or {"total_rows": 0})
    )

    metadata_payload: dict[str, Any] = {
        "table": resolved_table,
        "date_tag": date_tag,
        "parameters": dict(parameters or {}),
        "sources": _resolve_metadata_sources(pubchem_enabled),
        "outputs": outputs,
        "qc_summary": qc_payload,
    }
    if stats:
        metadata_payload["pipeline_stats"] = dict(stats)

    meta_path = write_meta_yaml(
        csv_path=dataset_path,
        stats=dict(stats) if stats is not None else None,
        extra_metadata=metadata_payload,
        generated_at=_metadata_generated_at(generated_at),
    )

    logger.info("[META] Метаданные сохранены: %s", meta_path)
    return meta_path


def _is_placeholder_user_agent(user_agent: str | None) -> bool:
    """Return ``True`` if *user_agent* still uses the default placeholder contact."""

    if not user_agent:
        return True
    return _PLACEHOLDER_CONTACT_EMAIL in user_agent


def _prepare_pubchem_api_cfg(cfg: Config, api_cfg: ApiCfg) -> ApiCfg:
    """Return an :class:`ApiCfg` with a valid PubChem user agent."""

    pubchem_user_agent = cfg.pubchem.user_agent.strip()
    api_user_agent = api_cfg.user_agent.strip()

    if not _is_placeholder_user_agent(pubchem_user_agent):
        if pubchem_user_agent != api_user_agent:
            return api_cfg.model_copy(update={"user_agent": pubchem_user_agent})
        return api_cfg

    if not _is_placeholder_user_agent(api_user_agent):
        return api_cfg

    raise ValueError(
        "PubChem configuration requires a user_agent with real contact details; "
        "set sources.pubchem.user_agent or api.user_agent in config.yaml.",
    )


def read_input_ids(
    input_csv: Path,
    *,
    column: str,
    io_cfg: IoCfg,
    limit: int | None,
    offset: int = 0,
) -> tuple[int, ReadInputIdsResult | None]:
    """Load identifiers from ``input_csv`` honouring ``limit`` when provided."""

    try:
        ids_iter = iter(
            io.read_ids(
                input_csv,
                column=column,
                cfg=io_cfg,
                keep_na_markers=io_cfg.keep_na_markers,
            )
        )
        if offset:
            ids_iter = islice(ids_iter, offset, None)
            logger.info("process_offset", offset=offset)
        if limit is not None:
            ids_iter = islice(ids_iter, limit)
        ids_iter, sample_iter = tee(ids_iter)
        sample_ids = tuple(islice(sample_iter, _FETCH_ERROR_SAMPLE_SIZE))
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "read_fail",
            error=str(exc),
            path=str(input_csv),
        )
        return 1, None

    return 0, ReadInputIdsResult(
        ids_iter=ids_iter,
        sample_ids=sample_ids,
    )


def _add_pipeline_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Attach acquisition metadata using the legacy pipeline helper."""

    from library.pipelines.common import add_pipeline_metadata as _add_pipeline_metadata

    return _add_pipeline_metadata(df)


def _log_missing_identifier_summary(identifiers: Sequence[str]) -> None:
    """Log a truncated list of missing identifiers including their total count."""

    if not identifiers:
        return

    sample = list(islice(identifiers, _MISSING_IDENTIFIER_LOG_SAMPLE_SIZE))
    logger.warning(
        "chembl_missing_identifiers",
        total=len(identifiers),
        sample=sample,
    )


def _emit_missing_identifier_logs(identifiers: Sequence[str]) -> None:
    """Emit warning and informational logs describing ``identifiers``."""

    if not identifiers:
        return

    _log_missing_identifier_summary(identifiers)
    logger.info(
        "chembl_missing_identifiers_total",
        count=len(identifiers),
    )


def generate_missing_identifier_placeholders(
    missing_ids: Sequence[str],
    *,
    columns: Sequence[str] | None = None,
    chunk_size: int = MISSING_IDENTIFIER_PLACEHOLDER_CHUNK_SIZE,
) -> Iterator[pd.DataFrame]:
    """Yield placeholder DataFrames for ``missing_ids`` in ``chunk_size`` batches."""

    if not missing_ids:
        return

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")

    if columns is None:
        columns = ("molecule_chembl_id",)

    ordered_columns = list(columns)
    if "molecule_chembl_id" not in ordered_columns:
        ordered_columns.append("molecule_chembl_id")

    total = len(missing_ids)
    for start in range(0, total, chunk_size):
        chunk = missing_ids[start : start + chunk_size]
        if not chunk:
            continue
        placeholder = pd.DataFrame(
            pd.NA,
            index=range(len(chunk)),
            columns=ordered_columns,
        )
        placeholder["molecule_chembl_id"] = chunk
        yield placeholder


def integrate_missing_identifiers(
    df: pd.DataFrame,
    *,
    missing_ids: Sequence[str],
    requested_ids: Sequence[str],
) -> pd.DataFrame:
    """Attach placeholder rows for ``missing_ids`` and restore input ordering."""

    if missing_ids:
        columns = list(df.columns)
        df = pd.concat(
            chain(
                [df],
                generate_missing_identifier_placeholders(
                    missing_ids, columns=columns
                ),
            ),
            ignore_index=True,
            sort=False,
        )
    else:
        df = df.reset_index(drop=True)

    if not requested_ids or "molecule_chembl_id" not in df.columns:
        return df

    order_map: dict[str, deque[int]] = {}
    for index, identifier in enumerate(requested_ids):
        if identifier is None:
            continue
        key = str(identifier).strip()
        if not key:
            continue
        upper_key = key.upper()
        order_map.setdefault(upper_key, deque()).append(index)

    if not order_map:
        return df

    fallback_index = len(requested_ids)

    def _pop_order(value: Any) -> int:
        if value is None or pd.isna(value):
            return fallback_index
        key = str(value).strip()
        if not key:
            return fallback_index
        upper_key = key.upper()
        queue = order_map.get(upper_key)
        if queue:
            return queue.popleft()
        return fallback_index

    order_series = df["molecule_chembl_id"].map(_pop_order)
    ordered = df.assign(_input_order=order_series).sort_values(
        "_input_order", kind="stable"
    )
    ordered = ordered.drop(columns="_input_order").reset_index(drop=True)
    return ordered


class TestitemFetchError(RuntimeError):
    """Raised when streaming retrieval of test item data fails."""


class TestitemPipelineStageError(RuntimeError):
    """Raised when a streaming pipeline stage returns a non-zero status."""

    def __init__(self, code: int, message: str | None = None) -> None:
        detail = message if message is not None else f"pipeline stage failed ({code})"
        super().__init__(detail)
        self.code = code


class StageExecutionBudget:
    """Track a shared execution budget for a pipeline stage."""

    def __init__(
        self,
        stage_name: str,
        *,
        minutes: float | None,
        logger: Any = logger,
    ) -> None:
        self.stage_name = stage_name
        self._logger = logger
        if minutes is None or minutes <= 0:
            self._budget_seconds: float | None = None
            self._deadline: float | None = None
        else:
            self._budget_seconds = float(minutes) * 60.0
            self._deadline = time.monotonic() + self._budget_seconds
        self._exhausted = False

    def start(self) -> None:
        if self._deadline is None or self._budget_seconds is None:
            return
        self._logger.info(
            f"{self.stage_name}_execution_budget_started",
            budget_seconds=int(self._budget_seconds),
        )

    def enforce(self, label: str | None = None) -> None:
        if self._deadline is None or self._budget_seconds is None:
            return
        remaining = self._deadline - time.monotonic()
        if remaining >= 0:
            return
        if not self._exhausted:
            self._logger.error(
                f"{self.stage_name}_execution_budget_exhausted",
                label=label,
                budget_seconds=int(self._budget_seconds),
            )
            self._exhausted = True
        raise TestitemPipelineStageError(
            1,
            (
                f"{self.stage_name} stage exceeded execution budget "
                f"({self._budget_seconds / 60:.1f} minutes)"
            ),
        )


class StageWatchdog:
    """Background timer that monitors progress for a pipeline stage."""

    def __init__(
        self,
        stage_name: str,
        *,
        idle_minutes: float,
        logger: Any = logger,
        check_interval: float = 60.0,
    ) -> None:
        self.stage_name = stage_name
        self._logger = logger
        self._idle_timeout_seconds = max(0.0, float(idle_minutes)) * 60.0
        self._check_interval = check_interval
        self._stop_event = threading.Event()
        self._timed_out = threading.Event()
        self._thread: threading.Thread | None = None
        self._last_activity = time.monotonic()
        self._effective_interval = 0.0

    def __enter__(self) -> StageWatchdog:
        self.start()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.stop()

    def start(self) -> None:
        if self._idle_timeout_seconds <= 0:
            return
        self._stop_event.clear()
        self._timed_out.clear()
        self._last_activity = time.monotonic()
        interval = self._idle_timeout_seconds / 2 if self._idle_timeout_seconds else 0
        self._effective_interval = max(
            1.0, min(self._check_interval, interval or self._check_interval)
        )
        self._thread = threading.Thread(
            target=self._monitor,
            name=f"{self.stage_name}-watchdog",
            daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        if self._thread is None:
            return
        self._stop_event.set()
        self._thread.join(timeout=1)
        self._thread = None

    def ping(self, event: str | None = None, **payload: Any) -> None:
        self._last_activity = time.monotonic()
        if event:
            self._logger.info(
                f"{self.stage_name}_watchdog_progress",
                watchdog_event=event,
                **payload,
            )

    def raise_if_timed_out(self) -> None:
        if self._timed_out.is_set():
            raise TestitemPipelineStageError(
                1,
                (
                    f"{self.stage_name} stage stalled (idle timeout "
                    f"{self._idle_timeout_seconds / 60:.1f} minutes)"
                ),
            )

    def _monitor(self) -> None:
        if self._idle_timeout_seconds <= 0:
            return
        interval = self._effective_interval or max(1.0, self._idle_timeout_seconds / 2)
        while not self._stop_event.wait(interval):
            idle = time.monotonic() - self._last_activity
            if idle >= self._idle_timeout_seconds:
                self._logger.error(
                    f"{self.stage_name}_watchdog_timeout",
                    idle_seconds=int(idle),
                    timeout_seconds=int(self._idle_timeout_seconds),
                )
                self._timed_out.set()
                return


@lru_cache(maxsize=1)
def _load_chembl_library():
    """Return the ChemBL integration module without creating import cycles."""

    from library.integration import chembl_library as chembl_lib

    return chembl_lib


@lru_cache(maxsize=1)
def _load_pipeline_metadata_adder():
    """Return the pipeline metadata helper lazily to avoid circular imports."""

    from library.pipelines.common import add_pipeline_metadata as add_metadata

    return add_metadata


@lru_cache(maxsize=1)
def _load_testitem_schema():
    """Return the schema model and normalizer lazily to avoid circular imports."""

    return TestitemsSchema, normalize_testitems


@lru_cache(maxsize=1)
def _load_pubchem_augmenter():
    """Return the PubChem augmentation helper lazily to avoid circular imports."""

    from .pubchem import augment_pubchem as _augment

    return _augment


def _batched(iterable: Iterable[str], size: int) -> Iterator[list[str]]:
    """Yield lists of at most ``size`` elements from ``iterable``."""

    chunk: list[str] = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _disabled_optional_columns(cfg: Config) -> frozenset[str]:
    """Return optional schema columns disabled by configuration."""

    disabled: set[str] = set()

    pubchem_cfg = getattr(cfg, "pubchem", None)
    pubchem_enabled = (
        True if pubchem_cfg is None else getattr(pubchem_cfg, "enable", True)
    )
    if not pubchem_enabled:
        disabled.update(_PUBCHEM_OPTIONAL_COLUMNS)

    enrichment_cfg = getattr(cfg, "testitem_molecule_enrichment", None)
    enrichment_enabled = (
        True if enrichment_cfg is None else getattr(enrichment_cfg, "enable", False)
    )
    if not enrichment_enabled:
        disabled.add(_SALT_OPTIONAL_COLUMN)

    return frozenset(disabled)


def fetch_testitems(
    ids_iter: Iterable[str],
    *,
    api_cfg: ApiCfg,
    batch_size: int,
    timeout: float,
    client: ChemblClient,
    sample_ids: Sequence[str],
    fields: Sequence[str] | None,
    page_limit: int,
    retry_cfg: TestitemBatchRetryCfg,
) -> tuple[int, Iterator[pd.DataFrame] | None, RequestedIdsSnapshot]:
    """Retrieve ChEMBL test item records for ``ids_iter`` in chunks."""

    logger.info("chembl_fetch_start", batch_size=batch_size)
    chembl_lib = _load_chembl_library()
    requested_ids_snapshot = RequestedIdsSnapshot()

    def _requested_iter() -> Iterator[str]:
        try:
            for identifier in ids_iter:
                requested_ids_snapshot.append(identifier)
                yield str(identifier)
        finally:
            requested_ids_snapshot.finish()
            if requested_ids_snapshot.duplicate_count:
                logger.warning(
                    "chembl_duplicate_requested_identifiers",
                    duplicate_count=requested_ids_snapshot.duplicate_count,
                    sample=list(requested_ids_snapshot.duplicate_sample),
                )

    seen_ids: set[str] = set()
    duplicate_ids: list[str] = []
    duplicate_ids_seen: set[str] = set()

    def _record_duplicates(values: Iterable[str]) -> None:
        for value in values:
            if value not in duplicate_ids_seen:
                duplicate_ids_seen.add(value)
                duplicate_ids.append(value)

    min_retry_size = max(1, retry_cfg.min_size) if retry_cfg.enable else 1

    def _should_retry(current_size: int) -> bool:
        if not retry_cfg.enable:
            return False
        return current_size > min_retry_size

    def _next_chunk_size(current_size: int) -> int:
        if not retry_cfg.enable:
            return current_size
        reduced = max(min_retry_size, int(current_size * retry_cfg.shrink_factor))
        if reduced >= current_size and current_size > min_retry_size:
            reduced = current_size - 1
        return max(min_retry_size, reduced)

    def _load_chunk(batch: Sequence[str], chunk_size: int) -> pd.DataFrame:
        try:
            frame = chembl_lib.get_testitem(
                batch,
                cfg=api_cfg,
                client=client,
                chunk_size=chunk_size,
                timeout=timeout,
                fields=fields,
                page_limit=page_limit,
            )
        except (
            requests.RequestException,
            ValueError,
        ) as exc:  # pragma: no cover - network
            if not _should_retry(chunk_size):
                logger.error(
                    "testitem_fetch_failed",
                    error=str(exc),
                    batch_size=chunk_size,
                    timeout=timeout,
                    sample_ids=list(sample_ids),
                )
                raise TestitemFetchError(str(exc)) from exc
            next_chunk_size = _next_chunk_size(chunk_size)
            if next_chunk_size >= chunk_size:
                logger.error(
                    "testitem_fetch_failed",
                    error=str(exc),
                    batch_size=chunk_size,
                    timeout=timeout,
                    sample_ids=list(sample_ids),
                )
                raise TestitemFetchError(str(exc)) from exc
            logger.warning(
                "testitem_fetch_retry_reduced_batch",
                error=str(exc),
                previous_chunk_size=chunk_size,
                next_chunk_size=next_chunk_size,
                identifier_count=len(batch),
            )
            return _load_chunk(batch, next_chunk_size)

        if fields:
            resolved_targets: dict[str, str] = {}
            for column in fields:
                if not isinstance(column, str):
                    continue
                target = TESTITEM_STRUCTURE_COLUMNS.get(column, column)
                resolved_targets[column] = target

            if resolved_targets:
                frame = frame.rename(columns=resolved_targets)

            missing_targets = [
                target
                for target in resolved_targets.values()
                if target not in frame.columns
            ]
            if missing_targets:
                na_series = pd.Series(pd.NA, index=frame.index, dtype="object")
                for target in missing_targets:
                    frame[target] = na_series.copy()

        if "molecule_chembl_id" not in frame.columns:
            frame["molecule_chembl_id"] = pd.Series(dtype="string")
        frame["molecule_chembl_id"] = (
            frame["molecule_chembl_id"].astype("string").str.upper()
        )

        duplicated_mask = frame["molecule_chembl_id"].duplicated(keep=False)
        if duplicated_mask.any():
            _record_duplicates(
                frame.loc[duplicated_mask, "molecule_chembl_id"]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

        cleaned = frame.where(pd.notna(frame), pd.NA).convert_dtypes()
        cleaned["molecule_chembl_id"] = cleaned["molecule_chembl_id"].astype("string")

        if seen_ids:
            keep_mask = ~cleaned["molecule_chembl_id"].isin(seen_ids)
            dropped = (
                cleaned.loc[~keep_mask, "molecule_chembl_id"]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            _record_duplicates(dropped)
            cleaned = cleaned.loc[keep_mask]

        if not cleaned.empty:
            seen_ids.update(
                cleaned["molecule_chembl_id"].dropna().astype(str).unique().tolist()
            )

        return cleaned.reset_index(drop=True)

    requested_iter = _requested_iter()
    batched_iter = _batched(requested_iter, max(1, batch_size))

    try:
        first_batch = next(batched_iter)
    except StopIteration:
        requested_iter.close()
        logger.info("chembl_fetch_done", rows=0)
        logger.info("identifiers_retrieved", count=0)
        if duplicate_ids:
            duplicates_truncated = (
                len(duplicate_ids) > _DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE
            )
            logger.warning(
                "chembl_duplicate_identifiers",
                duplicate_count=len(duplicate_ids),
                duplicate_ids=duplicate_ids[:_DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE],
                duplicates_truncated=duplicates_truncated,
            )
        return 0, iter(()), requested_ids_snapshot

    prefetched: list[pd.DataFrame] = []

    try:
        prefetched.append(_load_chunk(first_batch, batch_size))
    except TestitemFetchError:
        batched_iter.close()
        requested_iter.close()
        return 1, None, requested_ids_snapshot

    rows_counter = 0

    def _chunk_stream() -> Iterator[pd.DataFrame]:
        nonlocal rows_counter
        try:
            for prefetched_chunk in prefetched:
                if not prefetched_chunk.empty:
                    rows_counter += len(prefetched_chunk)
                    yield prefetched_chunk
            for batch in batched_iter:
                chunk_df = _load_chunk(batch, batch_size)
                if not chunk_df.empty:
                    rows_counter += len(chunk_df)
                    yield chunk_df
        finally:
            logger.info("chembl_fetch_done", rows=rows_counter)
            logger.info("identifiers_retrieved", count=rows_counter)
            if duplicate_ids:
                duplicates_truncated = (
                    len(duplicate_ids) > _DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE
                )
                logger.warning(
                    "chembl_duplicate_identifiers",
                    duplicate_count=len(duplicate_ids),
                    duplicate_ids=duplicate_ids[
                        : _DUPLICATE_IDENTIFIER_LOG_SAMPLE_SIZE
                    ],
                    duplicates_truncated=duplicates_truncated,
                )

    return 0, _chunk_stream(), requested_ids_snapshot


def apply_testitem_enrichment(
    df: pd.DataFrame,
    *,
    enrichment_cfg: TestitemMoleculeEnrichmentCfg,
    io_cfg: IoCfg,
) -> tuple[int, pd.DataFrame | None]:
    """Apply optional test item enrichment if enabled."""

    if not enrichment_cfg.enable:
        return 0, df

    logger.info("testitem_enrichment_start")
    try:
        enriched = testitem_enrichment.enrich(
            df,
            cfg=enrichment_cfg,
            io_cfg=io_cfg,
        )
    except ValueError as exc:
        logger.error("testitem_enrichment_failed", error=str(exc))
        return 1, None
    logger.info("testitem_enrichment_done")
    return 0, enriched


def run_testitem_pipeline(
    cfg: Config,
    options: TestitemPipelineOptions,
) -> tuple[int, io.StandardOutputArtifacts | None]:
    """Execute compound retrieval from the ChEMBL API and augment with PubChem data."""

    limit = options.limit if options.limit is not None else cfg.testitem.limit
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="testitem.limit", limit=limit)
        return 1

    api_overrides: dict[str, Any] = {}
    if cfg.testitem.retries is not None:
        api_overrides["retries"] = cfg.testitem.retries
    if cfg.testitem.backoff_factor is not None:
        api_overrides["backoff_factor"] = cfg.testitem.backoff_factor
    api_cfg = cfg.api.model_copy(update=api_overrides) if api_overrides else cfg.api

    pubchem_override = getattr(options, "pubchem_enabled", None)
    try:
        default_pubchem_enabled = getattr(cfg.pubchem, "enable", True)
    except AttributeError:
        logger.warning(
            "pubchem_configuration_missing",
            reason="missing_attribute",
            detail="cfg.pubchem.enable is not accessible; defaulting to enabled.",
        )
        default_pubchem_enabled = True

    if pubchem_override is not None:
        pubchem_enabled = bool(pubchem_override)
    else:
        pubchem_enabled = bool(default_pubchem_enabled)
    if not pubchem_enabled:
        logger.warning(
            "pubchem_augmentation_disabled",
            reason="config_disabled",
            detail="PubChem augmentation is disabled; pubchem_* columns will remain empty.",
        )
    else:
        logger.info("pubchem_augmentation_enabled")
    pubchem_api_cfg = api_cfg
    if pubchem_enabled:
        pubchem_api_cfg = _prepare_pubchem_api_cfg(cfg, api_cfg)
        pc.init_session(pubchem_api_cfg, cfg.retry)

    requested_ids: Sequence[str] = ()
    missing_ids: list[str] = []
    input_csv = Path(options.input_csv)
    output_csv = Path(options.output_csv) if options.output_csv is not None else None
    offset = options.offset if options.offset is not None else cfg.testitem.offset

    def _client_factory(context: ETLContext) -> ChemblClient:
        return ChemblClient(
            api_cfg,
            cfg.retry,
            cfg.chembl,
            global_limiter=context.global_limiter,
        )

    with ETLContext(cfg, chembl_client_factory=_client_factory) as context:
        client = context.chembl_client
        read_status, read_result = read_input_ids(
            input_csv,
            column=cfg.testitem.column,
            io_cfg=cfg.io,
            limit=limit,
            offset=offset,
        )
        if read_status != 0 or read_result is None:
            return read_status, None

        fetch_status, chunk_iter, requested_ids = fetch_testitems(
            read_result.ids_iter,
            api_cfg=api_cfg,
            batch_size=cfg.testitem.batch_size,
            timeout=cfg.testitem.timeout,
            client=client,
            sample_ids=read_result.sample_ids,
            fields=cfg.testitem.fields,
            page_limit=cfg.testitem.request_limit,
            retry_cfg=cfg.testitem.batch_retry,
        )
        if fetch_status != 0 or chunk_iter is None:
            return fetch_status, None

        requested_unique = requested_ids.unique_map

        fetched_ids: set[str] = set()

        missing_ids = []
        parent_stats_holder: dict[str, ParentLookupStats] = {
            "value": ParentLookupStats(
                source=PARENT_LOOKUP_SOURCE_SKIPPED,
                missing=0,
                unique=0,
                attached=0,
                uncovered=0,
            )
        }
        rows_counter = 0

        enrichment_sources = getattr(cfg.testitem_molecule_enrichment, "sources", None)
        hierarchy_lookup_path = (
            getattr(enrichment_sources, "molecule_hierarchy_path", None)
            if enrichment_sources is not None
            else None
        )
        parent_timeout = cfg.testitem.timeout
        pubchem_fallback_context = (
            PubChemAugmentationContext(
                pubchem_cfg=cfg.pubchem,
                api_cfg=pubchem_api_cfg,
                retry_cfg=cfg.retry,
                client=client,
                timeout=parent_timeout,
                fields=cfg.testitem.fields,
                request_limit=cfg.testitem.request_limit,
            )
            if pubchem_enabled
            else None
        )

        def _parent_stats_supplier() -> ParentLookupStats:
            return parent_stats_holder["value"]

        def _processed_chunks() -> Iterator[pd.DataFrame]:
            nonlocal rows_counter
            try:
                for chunk in chunk_iter:
                    rows_counter += len(chunk)
                    prep_status, prep = prepare_parent_enrichment(
                        chunk,
                        catalog_cfg=cfg.molecule_catalog,
                        io_cfg=cfg.io,
                        api_cfg=api_cfg,
                        timeout=parent_timeout,
                        client=client,
                        hierarchy_lookup_path=hierarchy_lookup_path,
                    )
                    if prep_status != 0 or prep is None:
                        raise TestitemPipelineStageError(prep_status)

                    parent_status, parent_result = run_parent_enrichment(
                        prep,
                        client=client,
                        api_cfg=api_cfg,
                        catalog_cfg=cfg.molecule_catalog,
                        timeout=parent_timeout,
                    )
                    if parent_status != 0 or parent_result is None:
                        raise TestitemPipelineStageError(parent_status)

                    current = parent_result.df
                    parent_stats_holder["value"] = _merge_parent_stats(
                        parent_stats_holder["value"], parent_result.parent_stats
                    )

                    if pubchem_enabled:
                        current = _load_pubchem_augmenter()(
                            current,
                            pubchem_cfg=cfg.pubchem,
                            api_cfg=pubchem_api_cfg,
                            retry_cfg=cfg.retry,
                            timeout=parent_timeout,
                            client=client,
                            fields=cfg.testitem.fields,
                            request_limit=cfg.testitem.request_limit,
                        )

                    enrichment_status, enriched_df = apply_testitem_enrichment(
                        current,
                        enrichment_cfg=cfg.testitem_molecule_enrichment,
                        io_cfg=cfg.io,
                    )
                    if enrichment_status != 0 or enriched_df is None:
                        raise TestitemPipelineStageError(enrichment_status)

                    if "molecule_chembl_id" in enriched_df.columns:
                        ids_series = (
                            enriched_df["molecule_chembl_id"]
                            .dropna()
                            .astype(str)
                            .str.strip()
                        )
                        fetched_ids.update(
                            {value.upper() for value in ids_series if value}
                        )

                    yield enriched_df
            except (
                TestitemFetchError
            ) as exc:  # pragma: no cover - propagated network error
                raise TestitemPipelineStageError(1, str(exc)) from exc

            for key, original in requested_unique.items():
                if key not in fetched_ids:
                    missing_ids.append(original)
            if missing_ids:
                _emit_missing_identifier_logs(missing_ids)
                for placeholder in generate_missing_identifier_placeholders(
                    missing_ids, columns=("molecule_chembl_id",)
                ):
                    rows_counter += len(placeholder)
                    yield placeholder

        output_path = (
            output_csv
            if output_csv is not None
            else io.default_output_path(
                input_csv,
                cfg.io,
                date=getattr(options, "date", None),
            )
        )

        artifacts: io.StandardOutputArtifacts | None = None
        try:
            exit_code, artifacts = finalize_output(
                _processed_chunks(),
                cfg=cfg,
                output=output_path,
                parent_stats_supplier=_parent_stats_supplier,
                input_csv=input_csv,
                missing_ids=missing_ids,
                emit_legacy_artifacts=options.emit_legacy_artifacts,
                pubchem_context=pubchem_fallback_context,
                options=options,
            )
        except TestitemPipelineStageError as exc:
            return exc.code, None

    if limit is not None:
        logger.info("process_limit", limit=min(limit, rows_counter))

    return exit_code, artifacts


def finalize_output(
    chunks: Iterable[pd.DataFrame],
    *,
    cfg: Config,
    output: Path,
    parent_stats_supplier: Callable[[], ParentLookupStats],
    input_csv: Path,
    missing_ids: Sequence[str] | None = None,
    emit_legacy_artifacts: bool = False,
    pubchem_context: PubChemAugmentationContext | None = None,
    options: TestitemPipelineOptions | None = None,
) -> tuple[int, io.StandardOutputArtifacts | None]:
    """Normalise, validate, and persist the final dataset from ``chunks``."""

    schema_model, normalizer = _load_testitem_schema()
    schema_cols = list(schema_model.columns)
    pubchem_cfg = getattr(cfg, "pubchem", None)
    pubchem_enabled = (
        True if pubchem_cfg is None else getattr(pubchem_cfg, "enable", True)
    )
    # ``pubchem_augmentation_enabled`` captures the configuration state so that
    # downstream telemetry remains stable even when no augmentation takes
    # place (for example when the dataset is empty).  The pipeline historically
    # emitted this flag under the ``pubchem_augmentation_enabled`` key, so keep
    # a dedicated variable instead of reusing ``pubchem_enabled`` directly.
    pubchem_augmentation_enabled = bool(pubchem_enabled)
    disabled_optional = _disabled_optional_columns(cfg)
    required_cols = {
        name
        for name, col in schema_model.columns.items()
        if col.required and name not in disabled_optional
    }
    optional_cols = (set(schema_model.columns) - required_cols) - disabled_optional
    col_order = schema_cols
    key_cols = ["molecule_chembl_id"]

    failure_chunk_size = resolve_failure_chunk_size(cfg)
    rows_total = 0
    rows_written = 0
    exit_code = 0
    next_raw_index = 0
    columns_seen: set[str] = set()
    columns_present: set[str] = set()
    columns_to_fill: set[str] = set(optional_cols)
    expected_columns: set[str] = set()
    column_dtypes: dict[str, pd.api.extensions.ExtensionDtype | str | type | None] = {}
    failure_cases = SidecarErrors(chunk_size=failure_chunk_size)
    failure_count = 0

    chunk_iter = iter(chunks)

    schema_dtype_lookup: dict[str, str] = {
        name: str(col.dtype) for name, col in schema_model.columns.items()
    }

    def _column_dtype(column: str) -> pd.api.extensions.ExtensionDtype | str | type:
        dtype_name = schema_dtype_lookup.get(column, "string")
        if dtype_name in {"str", "string"}:
            return pd.StringDtype()
        if dtype_name == "boolean":
            return pd.BooleanDtype()
        if dtype_name == "object":
            return object
        return pd.StringDtype()

    def _normalise_dtype(
        dtype: pd.api.extensions.ExtensionDtype | str | type | None,
    ) -> pd.api.extensions.ExtensionDtype | str | type:
        if dtype is None:
            return pd.StringDtype()
        try:
            pandas_dtype = pd.api.types.pandas_dtype(dtype)
        except TypeError:
            return dtype
        if isinstance(pandas_dtype, pd.BooleanDtype):
            return pd.BooleanDtype()
        if isinstance(pandas_dtype, pd.StringDtype):
            return pd.StringDtype()
        if pd.api.types.is_integer_dtype(pandas_dtype):
            return pd.Int64Dtype()
        if pd.api.types.is_float_dtype(pandas_dtype):
            return pd.Float64Dtype()
        return pandas_dtype

    def _ensure_column_alignment(frame: pd.DataFrame) -> None:
        missing = (expected_columns | columns_to_fill | columns_seen) - set(
            frame.columns
        )
        if not missing:
            return
        for column in sorted(missing):
            dtype = column_dtypes.get(column)
            if dtype is None:
                dtype = _column_dtype(column)
            dtype = _normalise_dtype(dtype)
            column_dtypes[column] = dtype
            if isinstance(dtype, pd.api.extensions.ExtensionDtype):
                empty_values = pd.array([pd.NA] * len(frame.index), dtype=dtype)
                frame[column] = pd.Series(empty_values, index=frame.index)
            else:
                placeholder = pd.Series(
                    [pd.NA] * len(frame.index), index=frame.index, dtype="object"
                )
                if dtype is not object:
                    placeholder = placeholder.astype(dtype, copy=False)
                frame[column] = placeholder

    def _process_chunk(raw: pd.DataFrame) -> pd.DataFrame:
        nonlocal rows_total, rows_written, exit_code, columns_seen, failure_count
        nonlocal next_raw_index

        rows_total += len(raw)
        current = normalizer(raw)
        if "pubchem_cid" in current.columns:
            current["pubchem_cid"] = current["pubchem_cid"].astype(object)
        current = _add_pipeline_metadata(current)
        ensure_raw_index_column(current)
        start_index = next_raw_index
        columns_present.update(current.columns)
        _ensure_column_alignment(current)
        for column in current.columns:
            column_dtypes.setdefault(column, _normalise_dtype(current.dtypes[column]))
        columns_seen.update(current.columns)
        expected_columns.update(columns_seen)
        _ensure_column_alignment(current)

        chunk_missing_required = required_cols - set(current.columns)
        if chunk_missing_required:
            logger.warning(
                "validation_skipped",
                missing_columns=sorted(chunk_missing_required),
            )
            raise TestitemPipelineStageError(
                1,
                ", ".join(sorted(chunk_missing_required)),
            )

        try:
            validation = validate_testitems(current, return_result=True)
        except SchemaErrors as exc:
            for row in exc.failure_cases.to_dict("records"):
                failure_cases.add_error(row)
                failure_count += 1
            exit_code = 1
            validated = getattr(exc, "validated_data", current)
        else:
            validated = validation.data
            if not validation.failure_cases.empty:
                exit_code = 1
                for row in validation.failure_cases.to_dict("records"):
                    failure_cases.add_error(row)
                    failure_count += 1

        end_index = start_index + len(validated)
        if RAW_INDEX_COLUMN in validated.columns:
            validated.loc[:, RAW_INDEX_COLUMN] = pd.Series(
                pd.RangeIndex(start_index, end_index),
                index=validated.index,
                dtype="int64",
                name=RAW_INDEX_COLUMN,
            )
        else:
            validated.insert(
                0,
                RAW_INDEX_COLUMN,
                pd.Series(
                    pd.RangeIndex(start_index, end_index),
                    index=validated.index,
                    dtype="int64",
                    name=RAW_INDEX_COLUMN,
                ),
            )
        next_raw_index = end_index

        rows_written += len(validated)
        return validated

    prepared_chunks: list[pd.DataFrame] = []

    try:
        first_raw = next(chunk_iter)
    except StopIteration:
        empty = pd.DataFrame(columns=["molecule_chembl_id"])
        try:
            prepared_chunks.append(_process_chunk(empty))
        except TestitemPipelineStageError:
            return 1, None
    else:
        try:
            prepared_chunks.append(_process_chunk(first_raw))
        except TestitemPipelineStageError:
            return 1, None

    def _validated_chunks() -> Iterator[pd.DataFrame]:
        for chunk in prepared_chunks:
            if not chunk.empty or columns_seen:
                _ensure_column_alignment(chunk)
                yield chunk
        for raw_chunk in chunk_iter:
            processed = _process_chunk(raw_chunk)
            _ensure_column_alignment(processed)
            yield processed

    missing_required = required_cols - columns_seen
    if missing_required:
        logger.warning(
            "validation_skipped",
            missing_columns=sorted(missing_required),
        )
        return 1, None

    if failure_count:
        failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")
        if emit_legacy_artifacts:
            failure_cases.save(failure_path, cfg=cfg)
            logger.error(
                "validation_failed",
                failures=failure_count,
                path=str(failure_path),
            )
        else:
            logger.error("validation_failed", failures=failure_count)

    if exit_code != 0:
        return exit_code, None

    pubchem_fallback_used = False
    pubchem_fallback_applied = False

    qc_ratio_counts: dict[str, int] = {column: 0 for column in _QC_RATIO_COLUMNS}
    pubchem_value_columns: set[str] = set()

    artifacts: io.StandardOutputArtifacts | None = None
    qc_summary: dict[str, Any] = {"total_rows": 0}

    with tempfile.TemporaryDirectory() as staging_dir:
        staging_dir_path = Path(staging_dir)
        chunk_paths: list[Path] = []
        fallback_candidates: list[pd.DataFrame] = []
        extra_columns_order: list[str] = []

        for index, chunk in enumerate(_validated_chunks()):
            if chunk.empty and not columns_seen:
                _ensure_column_alignment(chunk)
            for column in chunk.columns:
                if column not in col_order and column not in extra_columns_order:
                    extra_columns_order.append(column)
            if pubchem_context is not None and not chunk.empty:
                available_columns = [
                    column
                    for column in _PUBCHEM_OPTIONAL_COLUMNS
                    if column in chunk.columns
                ]
                if available_columns:
                    pubchem_columns = chunk[available_columns].replace("", pd.NA)
                    key_columns = [
                        column
                        for column in _PUBCHEM_KEY_COLUMNS
                        if column in pubchem_columns.columns
                    ]
                    if key_columns:
                        missing_mask = pubchem_columns[key_columns].isna().any(axis=1)
                    else:
                        missing_mask = pd.Series(True, index=chunk.index)
                else:
                    missing_mask = pd.Series(True, index=chunk.index)
                if missing_mask.any():
                    fallback_candidates.append(chunk.loc[missing_mask].copy())

            chunk_path = staging_dir_path / f"chunk_{index:06d}.parquet"
            chunk.to_parquet(chunk_path, index=False)
            chunk_paths.append(chunk_path)

        missing_optional = optional_cols - columns_present
        if missing_optional:
            logger.warning(
                "optional_columns_missing",
                columns=sorted(missing_optional),
            )

        raw_index_added = RAW_INDEX_COLUMN not in columns_present

        table_name, date_tag = io.derive_output_labels(
            output,
            default_table=_DEFAULT_TABLE_NAME,
            fallback_date=getattr(
                getattr(cfg, "io", None), "default_date_prefix", None
            ),
        )

        fallback_updates_frame: pd.DataFrame | None = None
        if pubchem_context is not None and fallback_candidates:
            logger.info("pubchem_fallback_augment_start")
            to_augment = pd.concat(fallback_candidates, ignore_index=False)
            augmented_subset = _load_pubchem_augmenter()(
                to_augment,
                pubchem_cfg=pubchem_context.pubchem_cfg,
                api_cfg=pubchem_context.api_cfg,
                retry_cfg=pubchem_context.retry_cfg,
                timeout=pubchem_context.timeout,
                client=pubchem_context.client,
                fields=pubchem_context.fields,
                request_limit=pubchem_context.request_limit,
            )
            if RAW_INDEX_COLUMN not in augmented_subset.columns:
                augmented_subset.insert(
                    0,
                    RAW_INDEX_COLUMN,
                    to_augment[RAW_INDEX_COLUMN].to_numpy(),
                )
            fallback_updates_frame = augmented_subset.set_index(RAW_INDEX_COLUMN)
            fallback_updates_frame.index = fallback_updates_frame.index.astype(
                "int64"
            )
            pubchem_fallback_used = True
            pubchem_fallback_applied = True
            logger.info("pubchem_fallback_augment_done")

        ordered_columns = [
            column for column in col_order if column in columns_seen
        ]
        extra_columns = [
            column for column in extra_columns_order if column not in ordered_columns
        ]
        export_column_order = ordered_columns + extra_columns

        logger.info(
            "raw_index_column_added",
            column=RAW_INDEX_COLUMN,
            rows=int(rows_written),
            inserted=raw_index_added,
        )

        doc_quality_cfg = cfg.system.doc_quality
        include_columns = getattr(doc_quality_cfg, "include_columns", None)
        exclude_columns = getattr(doc_quality_cfg, "exclude_columns", None)
        sample_rows_cfg = getattr(doc_quality_cfg, "sample_rows", None)
        auto_sample_limit = getattr(doc_quality_cfg, "auto_sample_row_limit", None)
        correlation_max_columns = getattr(
            doc_quality_cfg, "correlation_max_columns", None
        )

        row_count = int(rows_written)
        effective_sample_rows = sample_rows_cfg
        if (
            auto_sample_limit is not None
            and row_count > auto_sample_limit
            and (
                effective_sample_rows is None
                or effective_sample_rows > int(auto_sample_limit)
            )
        ):
            effective_sample_rows = int(auto_sample_limit)
            logger.warning(
                "doc_quality_sampling_applied",
                rows=row_count,
                sample_rows=effective_sample_rows,
                limit=int(auto_sample_limit),
            )

        correlation_include_columns = include_columns
        if (
            correlation_max_columns is not None
            and correlation_max_columns >= 1
            and include_columns is None
        ):
            exclude_set = set(exclude_columns or ())
            numeric_columns = [
                column
                for column in export_column_order
                if column not in exclude_set
                and column != RAW_INDEX_COLUMN
                and pd.api.types.is_numeric_dtype(column_dtypes.get(column))
            ]
            if len(numeric_columns) > int(correlation_max_columns):
                limited_columns = list(numeric_columns[: int(correlation_max_columns)])
                correlation_include_columns = tuple(limited_columns)
                logger.warning(
                    "correlation_columns_sampled",
                    total=len(numeric_columns),
                    limit=int(correlation_max_columns),
                    columns=limited_columns,
                )

        dataset_path = Path(output)
        dataset_path.parent.mkdir(parents=True, exist_ok=True)

        qc_profiler = qc_report.TableQualityProfiler()
        corr_profiler = qc_report.TableQualityProfiler()

        qc_include_tuple = (
            tuple(include_columns) if include_columns is not None else None
        )
        qc_exclude_tuple = (
            tuple(exclude_columns) if exclude_columns is not None else None
        )
        corr_include_tuple = (
            tuple(correlation_include_columns)
            if correlation_include_columns is not None
            else None
        )
        corr_exclude_tuple = qc_exclude_tuple

        qc_include_warning = [False]
        qc_exclude_warning = [False]
        qc_no_columns_warning = [False]
        corr_include_warning = [False]
        corr_exclude_warning = [False]
        corr_no_columns_warning = [False]

        remaining_qc_rows = effective_sample_rows
        remaining_corr_rows = effective_sample_rows

        def _apply_fallback(chunk: pd.DataFrame) -> pd.DataFrame:
            if fallback_updates_frame is None or fallback_updates_frame.empty:
                return chunk
            raw_indices = chunk[RAW_INDEX_COLUMN].astype("int64", copy=False)
            mask = raw_indices.isin(fallback_updates_frame.index)
            if not mask.any():
                return chunk
            replacements = fallback_updates_frame.reindex(raw_indices[mask])
            replacements.index = chunk.index[mask]
            updated = chunk.copy()
            for column in replacements.columns:
                updated.loc[mask, column] = replacements[column].to_numpy()
            return updated

        def _iter_final_chunks() -> Iterator[pd.DataFrame]:
            nonlocal remaining_qc_rows, remaining_corr_rows
            for path in chunk_paths:
                chunk = pd.read_parquet(path)
                chunk = _apply_fallback(chunk)
                chunk = chunk.reindex(columns=export_column_order, copy=False)

                for column in _QC_RATIO_COLUMNS:
                    if column in chunk.columns:
                        qc_ratio_counts[column] += int(chunk[column].notna().sum())

                for column in _PUBCHEM_OPTIONAL_COLUMNS:
                    if column in chunk.columns and column not in pubchem_value_columns:
                        series = chunk[column]
                        if not series.replace("", pd.NA).isna().all():
                            pubchem_value_columns.add(column)

                filtered_qc, remaining_qc_rows = _apply_sampling_and_filters(
                    chunk,
                    table_name=table_name,
                    sample_rows=remaining_qc_rows,
                    include_columns=qc_include_tuple,
                    exclude_columns=qc_exclude_tuple,
                    include_warning_logged=qc_include_warning,
                    exclude_warning_logged=qc_exclude_warning,
                    no_columns_logged=qc_no_columns_warning,
                )
                qc_profiler.consume(filtered_qc)

                filtered_corr, remaining_corr_rows = _apply_sampling_and_filters(
                    chunk,
                    table_name=table_name,
                    sample_rows=remaining_corr_rows,
                    include_columns=corr_include_tuple,
                    exclude_columns=corr_exclude_tuple,
                    include_warning_logged=corr_include_warning,
                    exclude_warning_logged=corr_exclude_warning,
                    no_columns_logged=corr_no_columns_warning,
                )
                corr_profiler.consume(filtered_corr)

                yield chunk

        try:
            dataset_path = io.write_csv(
                _iter_final_chunks(),
                dataset_path,
                cfg=cfg,
                key_cols=key_cols,
                col_order=export_column_order,
            )

            for sidecar in (
                dataset_path.with_suffix(dataset_path.suffix + ".meta.yaml"),
                dataset_path.with_suffix(dataset_path.suffix + ".meta.yaml.lock"),
            ):
                sidecar.unlink(missing_ok=True)

            quality_report = qc_report.build_qc_summary(
                None,
                table_name=table_name,
                include_columns=include_columns,
                exclude_columns=exclude_columns,
                sample_rows=effective_sample_rows,
                profiler=qc_profiler,
            )
            correlation_report = data_correlation.build_correlation_matrix(
                None,
                table_name=table_name,
                include_columns=correlation_include_columns,
                exclude_columns=exclude_columns,
                sample_rows=effective_sample_rows,
                profiler=corr_profiler,
            )

            correlation_path = dataset_path.parent / (
                f"{dataset_path.stem}_data_correlation_report_table.csv"
            )
            quality_path = dataset_path.parent / (
                f"{dataset_path.stem}_quality_report_table.csv"
            )

            corr_key_cols: list[str] = []
            if not correlation_report.empty:
                corr_key_cols = [str(correlation_report.columns[0])]
            io.write_csv(
                correlation_report,
                correlation_path,
                cfg=cfg,
                key_cols=corr_key_cols,
                col_order=list(correlation_report.columns),
            )

            qc_key_cols: list[str] = []
            if not quality_report.empty:
                qc_key_cols = [str(quality_report.columns[0])]
            io.write_csv(
                quality_report,
                quality_path,
                cfg=cfg,
                key_cols=qc_key_cols,
                col_order=list(quality_report.columns),
            )

            for sidecar in (
                correlation_path.with_suffix(correlation_path.suffix + ".meta.yaml"),
                correlation_path.with_suffix(
                    correlation_path.suffix + ".meta.yaml.lock"
                ),
                quality_path.with_suffix(quality_path.suffix + ".meta.yaml"),
                quality_path.with_suffix(quality_path.suffix + ".meta.yaml.lock"),
            ):
                sidecar.unlink(missing_ok=True)

            artifacts = io.StandardOutputArtifacts(
                dataset=dataset_path,
                correlation_report=correlation_path,
                quality_report=quality_path,
            )
        except (OSError, ValueError) as exc:
            logger.error("write_fail", error=str(exc), path=str(output))
            return 1, None

        qc_summary = {"total_rows": row_count}
        if row_count > 0:
            ratios: dict[str, float] = {}
            for column in _QC_RATIO_COLUMNS:
                if column in export_column_order:
                    ratios[column] = round(
                        float(qc_ratio_counts.get(column, 0)) / float(row_count),
                        4,
                    )
            if ratios:
                qc_summary["non_null_ratio"] = ratios

    logger.info(
        "write_done",
        rows=rows_written,
        dataset=str(artifacts.dataset) if artifacts else str(output),
        quality_report=str(artifacts.quality_report) if artifacts else "",
        correlation_report=str(artifacts.correlation_report) if artifacts else "",
    )

    if not emit_legacy_artifacts:
        for path in (
            artifacts.dataset,
            artifacts.quality_report,
            artifacts.correlation_report,
        ):
            Path(f"{path}.meta.yaml").unlink(missing_ok=True)

    rows_dropped = rows_total - rows_written
    parent_stats = parent_stats_supplier()
    missing_ids_tuple = tuple(missing_ids or ())

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_written,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(artifacts.dataset),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
        "pubchem_augmentation_enabled": pubchem_augmentation_enabled,
        "parent_lookup_hierarchy_attached": parent_stats.hierarchy_attached,
        "parent_lookup_fallback_attached": parent_stats.fallback_attached,
        "parent_lookup_no_parent": parent_stats.no_parent,
    }

    pubchem_columns_present = [
        column
        for column in _PUBCHEM_OPTIONAL_COLUMNS
        if column in export_column_order
    ]
    pubchem_columns_present.sort()
    pubchem_columns_with_values = sorted(pubchem_value_columns)
    if (
        pubchem_enabled
        and rows_written
        and pubchem_columns_present
        and not pubchem_columns_with_values
    ):
        logger.warning(
            "pubchem_augmentation_missing_values",
            columns=pubchem_columns_present,
        )

    stats.update(
        {
            "pubchem_enabled": pubchem_enabled,
            "pubchem_columns_present": pubchem_columns_present,
            "pubchem_columns_with_values": pubchem_columns_with_values,
            "pubchem_values_present": bool(pubchem_columns_with_values),
            "pubchem_fallback_applied": pubchem_fallback_applied,
            "pubchem_fallback_used": pubchem_fallback_used,
        }
    )

    if missing_ids_tuple:
        _log_missing_identifier_summary(missing_ids_tuple)
        missing_ids_sample = list(
            islice(missing_ids_tuple, _MISSING_IDENTIFIER_STATS_SAMPLE_SIZE)
        )
        total_missing = len(missing_ids_tuple)
        stats["missing_molecule_ids"] = missing_ids_sample
        stats["missing_ids_sample"] = missing_ids_sample
        stats["missing_molecule_ids_count"] = total_missing
        stats["missing_molecule_ids_total"] = total_missing
        stats["missing_molecule_ids_truncated"] = len(missing_ids_sample) < total_missing
    if parent_stats.failed_ids:
        stats["parent_lookup_failed_ids"] = list(parent_stats.failed_ids)
        stats["parent_lookup_failed_count"] = len(parent_stats.failed_ids)
        logger.info(
            "parent_lookup_failures_recorded",
            count=len(parent_stats.failed_ids),
            meta_key="parent_lookup_failed_ids",
        )

    logger.info("testitem_stats", **stats)

    parameters = _extract_metadata_parameters(
        cfg,
        options,
        emit_legacy_artifacts=emit_legacy_artifacts,
        pubchem_enabled=pubchem_augmentation_enabled,
    )

    legacy_meta_path: Path | None = None
    if emit_legacy_artifacts:
        legacy_meta_temp = write_meta_yaml(
            csv_path=artifacts.dataset,
            command=" ".join(sys.argv),
            config=_serialize_paths(cfg.to_dict()),
            inputs={"input_csv": str(input_csv)},
            stats=stats,
            schema="TestitemsSchema",
        )
        desired_legacy_path = legacy_meta_temp.with_name(
            legacy_meta_temp.name.replace(".meta.yaml", ".legacy.meta.yaml")
        )
        try:
            legacy_content = legacy_meta_temp.read_text(encoding="utf-8")
            with desired_legacy_path.open("w", encoding="utf-8") as handle:
                handle.write(legacy_content)
        except OSError as exc:  # pragma: no cover - defensive guard
            logger.warning(
                "testitem_legacy_metadata_copy_failed",
                error=str(exc),
                source=str(legacy_meta_temp),
                destination=str(desired_legacy_path),
            )
            legacy_meta_path = legacy_meta_temp
        else:
            legacy_meta_path = desired_legacy_path
        # ``_write_primary_metadata`` overwrites the canonical path next.

    primary_meta_path = _write_primary_metadata(
        dataset_frame=None,
        dataset_path=artifacts.dataset,
        quality_path=artifacts.quality_report,
        correlation_path=artifacts.correlation_report,
        table_name=table_name,
        date_tag=date_tag,
        parameters=parameters,
        stats=stats,
        pubchem_enabled=pubchem_augmentation_enabled,
        qc_summary=qc_summary,
        generated_at=_metadata_generated_at(),
    )

    if not emit_legacy_artifacts:
        return exit_code, artifacts

    fatal_quality_error = bool(getattr(doc_quality_cfg, "fatal_on_error", False))
    table_quality = build_table_quality_hook(
        doc_quality_cfg,
        table_name=artifacts.dataset.with_suffix(""),
        destination=artifacts.dataset.parent,
    )
    try:
        table_quality(artifacts.dataset)
    except Exception as exc:
        tb = traceback.format_exc()
        record_quality_failure(
            primary_meta_path,
            error=str(exc),
            error_type=exc.__class__.__name__,
            traceback=tb,
            fatal=fatal_quality_error,
        )
        if legacy_meta_path and legacy_meta_path != primary_meta_path:
            record_quality_failure(
                legacy_meta_path,
                error=str(exc),
                error_type=exc.__class__.__name__,
                traceback=tb,
                fatal=fatal_quality_error,
            )
        log_kwargs = {
            "error": str(exc),
            "error_type": exc.__class__.__name__,
            "path": str(artifacts.dataset),
            "traceback": tb,
        }
        if fatal_quality_error:
            log_kwargs["fatal"] = True
        logger.warning("quality_report_failed", **log_kwargs)
        if fatal_quality_error:
            return 1, None

    return exit_code, artifacts
