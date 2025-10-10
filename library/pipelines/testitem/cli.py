"""Reusable components for the test item data acquisition pipeline."""

from __future__ import annotations

import sys
import threading
import time
import traceback
from collections import OrderedDict, deque
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import chain, islice
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import io
from library.clients import pubchem as pc
from library.common.csv_utils import write_csv_chunks_deterministic
from library.common.log import logger
from library.common.metadata import (
    Stats,
    file_sha256,
    record_quality_failure,
    write_meta_yaml,
)
from library.common.sidecar import SidecarErrors
from library.config import (
    ApiCfg,
    Config,
    IoCfg,
    TestitemBatchRetryCfg,
    TestitemMoleculeEnrichmentCfg,
    _serialize_paths,
)
from library.integration.chembl_client import ChemblClient
from library.orchestration import ETLContext
from library.qa.reporting import build_table_quality_hook
from library.qa.validation import validate_testitems
from library.schemas import TestitemsSchema, normalize_testitems

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
_MISSING_IDENTIFIER_LOG_SAMPLE_SIZE = 10
_PLACEHOLDER_CONTACT_EMAIL = "contact@example.org"


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
        ids_iter = io.read_ids(
            input_csv,
            column=column,
            cfg=io_cfg,
            keep_na_markers=io_cfg.keep_na_markers,
        )
        ids_iter = iter(ids_iter)
        if offset:
            ids_iter = islice(ids_iter, offset, None)
            ids_iter = iter(ids_iter)
            logger.info("process_offset", offset=offset)
        if limit is not None:
            ids_iter = islice(ids_iter, limit)
            ids_iter = iter(ids_iter)
        sample_ids = tuple(islice(ids_iter, _FETCH_ERROR_SAMPLE_SIZE))
        ids_iter = chain(sample_ids, ids_iter)
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


def integrate_missing_identifiers(
    df: pd.DataFrame,
    *,
    missing_ids: Sequence[str],
    requested_ids: Sequence[str],
) -> pd.DataFrame:
    """Attach placeholder rows for ``missing_ids`` and restore input ordering."""

    if missing_ids:
        columns = list(df.columns)
        if "molecule_chembl_id" not in columns:
            columns.append("molecule_chembl_id")
        missing_frame = pd.DataFrame(
            pd.NA,
            index=range(len(missing_ids)),
            columns=columns,
        )
        missing_frame["molecule_chembl_id"] = list(missing_ids)
        df = pd.concat([df, missing_frame], ignore_index=True, sort=False)
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
) -> tuple[int, Iterator[pd.DataFrame] | None, tuple[str, ...]]:
    """Retrieve ChEMBL test item records for ``ids_iter`` in chunks."""

    logger.info("chembl_fetch_start", batch_size=batch_size)
    chembl_lib = _load_chembl_library()
    requested_ids_original: list[str] = []

    for identifier in ids_iter:
        text_identifier = str(identifier)
        requested_ids_original.append(text_identifier)

    requested_iter = iter(requested_ids_original)

    seen_ids: set[str] = set()
    duplicate_ids: set[str] = set()

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
            missing_fields = [
                column
                for column in fields
                if isinstance(column, str) and column not in frame.columns
            ]
            if missing_fields:
                na_series = pd.Series(pd.NA, index=frame.index, dtype="object")
                for column in missing_fields:
                    frame[column] = na_series.copy()

        if "molecule_chembl_id" not in frame.columns:
            frame["molecule_chembl_id"] = pd.Series(dtype="string")
        frame["molecule_chembl_id"] = (
            frame["molecule_chembl_id"].astype("string").str.upper()
        )

        duplicated_mask = frame["molecule_chembl_id"].duplicated(keep=False)
        if duplicated_mask.any():
            duplicate_ids.update(
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
            duplicate_ids.update(dropped)
            cleaned = cleaned.loc[keep_mask]

        if not cleaned.empty:
            seen_ids.update(
                cleaned["molecule_chembl_id"].dropna().astype(str).unique().tolist()
            )

        return cleaned.reset_index(drop=True)

    batched_iter = _batched(requested_iter, max(1, batch_size))

    try:
        first_batch = next(batched_iter)
    except StopIteration:
        logger.info("chembl_fetch_done", rows=0)
        logger.info("identifiers_retrieved", count=0)
        if duplicate_ids:
            logger.warning(
                "chembl_duplicate_identifiers",
                duplicate_count=len(duplicate_ids),
                duplicate_ids=sorted(duplicate_ids),
            )
        return 0, iter(()), tuple(requested_ids_original)

    prefetched: list[pd.DataFrame] = []

    try:
        prefetched.append(_load_chunk(first_batch, batch_size))
    except TestitemFetchError:
        return 1, None, tuple()

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
                logger.warning(
                    "chembl_duplicate_identifiers",
                    duplicate_count=len(duplicate_ids),
                    duplicate_ids=sorted(duplicate_ids),
                )

    return 0, _chunk_stream(), tuple(requested_ids_original)


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
) -> int:
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

    pubchem_enabled = getattr(cfg.pubchem, "enable", True)
    pubchem_api_cfg = api_cfg
    if pubchem_enabled:
        pubchem_api_cfg = _prepare_pubchem_api_cfg(cfg, api_cfg)
        pc.init_session(pubchem_api_cfg, cfg.retry)

    requested_ids: tuple[str, ...] = ()
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
            return read_status

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
            return fetch_status

        requested_unique: OrderedDict[str, str] = OrderedDict()
        for identifier in requested_ids:
            key = identifier.strip()
            if not key:
                continue
            upper_key = key.upper()
            requested_unique.setdefault(upper_key, identifier)

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

            missing_keys = [key for key in requested_unique if key not in fetched_ids]
            missing_ids.extend(requested_unique[key] for key in missing_keys)
            if missing_ids:
                logger.warning(
                    "chembl_missing_identifiers",
                    count=len(missing_ids),
                    identifiers=missing_ids,
                )
                placeholder = pd.DataFrame({"molecule_chembl_id": list(missing_ids)})
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

        try:
            exit_code = finalize_output(
                _processed_chunks(),
                cfg=cfg,
                output=output_path,
                parent_stats_supplier=_parent_stats_supplier,
                input_csv=input_csv,
                missing_ids=missing_ids,
            )
        except TestitemPipelineStageError as exc:
            return exc.code

    if limit is not None:
        logger.info("process_limit", limit=min(limit, rows_counter))

    return exit_code


def finalize_output(
    chunks: Iterable[pd.DataFrame],
    *,
    cfg: Config,
    output: Path,
    parent_stats_supplier: Callable[[], ParentLookupStats],
    input_csv: Path,
    missing_ids: Sequence[str] | None = None,
) -> int:
    """Normalise, validate, and persist the final dataset from ``chunks``."""

    schema_model, normalizer = _load_testitem_schema()
    schema_cols = list(schema_model.columns)
    required_cols = {name for name, col in schema_model.columns.items() if col.required}
    optional_cols = set(schema_model.columns) - required_cols
    col_order = schema_cols
    key_cols = ["molecule_chembl_id"]
    csv_chunksize = cfg.io.csv_chunksize

    rows_total = 0
    rows_written = 0
    exit_code = 0
    columns_seen: set[str] = set()
    columns_present: set[str] = set()
    columns_to_fill: set[str] = set()
    expected_columns: set[str] = set()
    column_dtypes: dict[str, pd.api.extensions.ExtensionDtype | str | type | None] = {}
    failure_cases = SidecarErrors()
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

        rows_total += len(raw)
        current = normalizer(raw)
        if "pubchem_cid" in current.columns:
            current["pubchem_cid"] = current["pubchem_cid"].astype(object)
        current = _add_pipeline_metadata(current)
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
            return 1
    else:
        try:
            prepared_chunks.append(_process_chunk(first_raw))
        except TestitemPipelineStageError:
            return 1

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
        return 1

    validated_chunks_list = list(_validated_chunks())

    missing_optional = optional_cols - columns_present
    if missing_optional:
        columns_to_fill.update(missing_optional)
        expected_columns.update(columns_to_fill)
        for column in missing_optional:
            column_dtypes.setdefault(column, _column_dtype(column))
        for chunk in validated_chunks_list:
            _ensure_column_alignment(chunk)
        columns_seen.update(columns_to_fill)
        logger.warning(
            "optional_columns_missing",
            columns=sorted(missing_optional),
        )

    if failure_count:
        failure_path = Path(output).with_name(f"{Path(output).stem}_failure_cases.csv")
        failure_cases.save(failure_path, cfg=cfg)
        logger.error(
            "validation_failed",
            failures=failure_count,
            path=str(failure_path),
        )

    if exit_code != 0:
        return exit_code

    try:
        csv_path = write_csv_chunks_deterministic(
            validated_chunks_list,
            output,
            col_order=col_order,
            key_cols=key_cols,
            chunksize=csv_chunksize,
            sort_chunksize=csv_chunksize,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
        )
        logger.info("write_done", rows=rows_written, path=str(csv_path))
    except TestitemPipelineStageError:
        return 1
    except (OSError, ValueError) as exc:
        logger.error("write_fail", error=str(exc), path=str(output))
        return 1

    rows_dropped = rows_total - rows_written
    parent_stats = parent_stats_supplier()
    missing_ids_tuple = tuple(missing_ids or ())

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_written,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
        "parent_lookup_source": parent_stats.source,
        "parent_lookup_missing": parent_stats.missing,
        "parent_lookup_hierarchy_attached": parent_stats.hierarchy_attached,
        "parent_lookup_fallback_attached": parent_stats.fallback_attached,
        "parent_lookup_no_parent": parent_stats.no_parent,
    }
    if missing_ids_tuple:
        stats["missing_molecule_ids"] = list(missing_ids_tuple)
        stats["missing_molecule_ids_count"] = len(missing_ids_tuple)
    if parent_stats.failed_ids:
        stats["parent_lookup_failed_ids"] = list(parent_stats.failed_ids)
        stats["parent_lookup_failed_count"] = len(parent_stats.failed_ids)
        logger.info(
            "parent_lookup_failures_recorded",
            count=len(parent_stats.failed_ids),
            meta_key="parent_lookup_failed_ids",
        )

    meta_path = write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(input_csv)},
        stats=stats,
        schema="TestitemsSchema",
    )

    doc_quality_cfg = cfg.system.doc_quality
    fatal_quality_error = bool(getattr(doc_quality_cfg, "fatal_on_error", False))
    table_quality = build_table_quality_hook(
        doc_quality_cfg,
        table_name=output.with_suffix(""),
        destination=output.parent,
    )
    try:
        table_quality(csv_path)
    except Exception as exc:
        tb = traceback.format_exc()
        record_quality_failure(
            meta_path,
            error=str(exc),
            error_type=exc.__class__.__name__,
            traceback=tb,
            fatal=fatal_quality_error,
        )
        log_kwargs = {
            "error": str(exc),
            "error_type": exc.__class__.__name__,
            "path": str(output),
            "traceback": tb,
        }
        if fatal_quality_error:
            log_kwargs["fatal"] = True
        logger.warning("quality_report_failed", **log_kwargs)
        if fatal_quality_error:
            return 1

    return exit_code
