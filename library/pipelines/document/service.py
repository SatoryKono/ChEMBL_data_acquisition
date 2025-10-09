"""Service layer for orchestrating document pipeline data fetching."""

from __future__ import annotations

import difflib
import heapq
import weakref

from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from contextlib import AbstractContextManager, ExitStack, contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from threading import Lock, local
from typing import Any, Hashable, TypeVar

import pandas as pd
import requests

from library.common.log import logger
from library.common.rate_limiter import RateLimiter, get_global_limiter, get_limiter
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
    session_with_retry,
)
from library.integration import openalex_crossref_library as ocl
from library.integration import pubmed_library as pl
from library.integration import semantic_scholar_library as ssl
from library.pipelines.document.pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    merge_metadata,
    normalise_doi,
)
from library.clients import _chunked


T = TypeVar("T")


@dataclass
class FallbackDoiMetrics:
    """Track statistics collected while applying fallback DOI mappings."""

    csv_records: int = 0
    total_candidates: int = 0
    matched_pmids: set[str] = field(default_factory=set)
    applied: int = 0
    conflicts: int = 0
    samples: list[dict[str, object]] = field(default_factory=list)

    def mark_csv_records(self, count: int) -> None:
        self.csv_records = max(count, 0)

    def mark_total_candidates(self, count: int) -> None:
        self.total_candidates = max(count, 0)

    def record_match(self, pmid: str) -> None:
        self.matched_pmids.add(pmid)

    def record_conflict(self) -> None:
        self.conflicts += 1

    def record_application(self, sample: dict[str, object]) -> None:
        self.applied += 1
        if len(self.samples) < 5:
            self.samples.append(sample)

    def as_log_kwargs(self) -> dict[str, object]:
        matched = len(self.matched_pmids)
        missing = max(self.total_candidates - matched, 0)
        return {
            "csv_records": self.csv_records,
            "matched_pmids": matched,
            "applied": self.applied,
            "conflicts": self.conflicts,
            "missing": missing,
            "samples": self.samples,
        }


@dataclass
class FallbackDoiState:
    """Hold runtime information about fallback DOI configuration."""

    path: Path
    mapping: dict[str, str]
    metrics: FallbackDoiMetrics
    overwrite: bool


class DocumentPipeline:
    """Encapsulate document pipeline orchestration and helper routines."""

    DOCUMENT_PROGRESS_INFO_INTERVAL = 100
    DOCUMENT_FRAME_CONCAT_STRIDE = 16

    def __init__(self, cfg: Config | None = None) -> None:
        self.cfg = cfg or Config()

    @staticmethod
    def limit_iterable(
        iterable: Iterable[T],
        limit: int,
    ) -> tuple[Iterator[T], Callable[[], int]]:
        """Return an iterator capped at ``limit`` items and a counter callback."""

        if limit < 0:
            msg = "limit must be non-negative"
            raise ValueError(msg)

        source = iter(iterable)
        count = 0

        def _generator() -> Iterator[T]:
            nonlocal count
            for item in source:
                if count >= limit:
                    break
                count += 1
                yield item

        limited_iter = _generator()

        def _get_count() -> int:
            return count

        return limited_iter, _get_count

    @staticmethod
    def read_csv_chunks(
        path: Path,
        *,
        cfg: Config,
        chunk_size: int,
    ) -> Iterator[pd.DataFrame]:
        """Yield DataFrame chunks from ``path`` while ensuring the reader closes."""

        reader = pd.read_csv(
            path,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            chunksize=chunk_size,
        )
        try:
            yield from reader
        finally:
            reader.close()

    @staticmethod
    def build_missing_input_context(path: Path) -> dict[str, object]:
        """Return structured hints describing how to resolve a missing input file."""

        info: dict[str, object] = {}
        parent = path.parent if path.parent != Path("") else Path(".")
        info["hint"] = (
            "Populate the expected input file or update --input to point to an "
            "existing CSV"
        )

        if not parent.exists():
            info["missing_parent"] = str(parent)
            return info

        try:
            entries = sorted(p for p in parent.iterdir() if p.is_file())
        except OSError:
            entries = []

        if not entries:
            info["available_files"] = []
            return info

        names = [p.name for p in entries]
        matches = difflib.get_close_matches(path.name, names, n=5, cutoff=0.55)
        if matches:
            suggestions = [parent / name for name in matches]
            info["suggestions"] = [str(path) for path in suggestions]
            best_match = suggestions[0]
            info["did_you_mean"] = str(best_match)
            info["cli_hint"] = f"--input \"{best_match}\""
            return info

        same_suffix = [str(p) for p in entries if p.suffix == path.suffix and p.suffix]
        if same_suffix:
            info["available_with_same_suffix"] = same_suffix[:5]
            return info

        info["available_files"] = [str(p) for p in entries[:5]]
        return info

    @staticmethod
    def build_fallback_doi_map(
        frame: pd.DataFrame,
        *,
        pmid_column: str,
        doi_column: str,
        metrics: FallbackDoiMetrics | None = None,
    ) -> dict[str, str]:
        """Return a mapping of PubMed identifiers to DOI overrides."""

        if pmid_column not in frame.columns or doi_column not in frame.columns:
            missing = [
                column
                for column in (pmid_column, doi_column)
                if column not in frame.columns
            ]
            raise ValueError(
                f"missing columns in fallback DOI file: {', '.join(missing)}"
            )

        if metrics is not None:
            metrics.mark_csv_records(len(frame.index))

        mapping: dict[str, str] = {}
        for pmid_value, doi_value in frame[[pmid_column, doi_column]].itertuples(
            index=False
        ):
            if pd.isna(pmid_value):
                pmid = ""
            else:
                pmid = str(pmid_value).strip()
            if pd.isna(doi_value):
                doi = ""
            else:
                doi = normalise_doi(doi_value)
            if not pmid or not doi:
                continue
            mapping[pmid] = doi
        if metrics is not None:
            metrics.mark_total_candidates(len(mapping))
        return mapping

    @staticmethod
    def _collect_doi_sources(frame: pd.DataFrame, index: Hashable) -> list[str]:
        """Return human readable labels describing DOI origins for ``frame`` row."""

        columns_to_sources = (
            ("PubMed.DOI", "pubmed"),
            ("scholar.DOI", "semantic_scholar"),
            ("OpenAlex.DOI", "openalex"),
            ("crossref.DOI", "crossref"),
            ("ChEMBL.doi", "chembl"),
            ("doi", "canonical"),
        )
        sources: list[str] = []
        seen: set[str] = set()
        for column, label in columns_to_sources:
            if column not in frame.columns:
                continue
            value = frame.at[index, column]
            if not normalise_doi(value):
                continue
            if label in seen:
                continue
            sources.append(label)
            seen.add(label)
        return sources

    @classmethod
    def apply_fallback_doi(
        cls,
        frame: pd.DataFrame,
        *,
        fallback_map: Mapping[str, str] | None,
        overwrite: bool,
        metrics: FallbackDoiMetrics | None = None,
        pmid_column: str = "PubMed.PMID",
    ) -> pd.DataFrame:
        """Apply fallback DOI mapping to ``frame`` and update ``metrics``."""

        if not fallback_map or pmid_column not in frame.columns:
            return frame

        pmid_series = (
            frame[pmid_column]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        fallback_series = pmid_series.map(fallback_map).fillna("")
        if not (fallback_series != "").any():
            return frame

        canonical_columns = [col for col in ("doi", "ChEMBL.doi") if col in frame.columns]
        if not canonical_columns:
            return frame

        for idx, fallback_doi in fallback_series.items():
            if not fallback_doi:
                continue
            pmid = pmid_series.at[idx]
            if not pmid:
                continue
            if metrics is not None:
                metrics.record_match(pmid)
            sources = cls._collect_doi_sources(frame, idx)
            existing_doi = ""
            for column in canonical_columns:
                current = normalise_doi(frame.at[idx, column])
                if current:
                    existing_doi = current
                    break
            conflict = bool(existing_doi and existing_doi != fallback_doi)
            if conflict and metrics is not None:
                metrics.record_conflict()
            if existing_doi and not overwrite:
                continue
            if existing_doi == fallback_doi:
                continue
            for column in canonical_columns:
                frame.at[idx, column] = fallback_doi
            if "doi_normalised" in frame.columns:
                frame.at[idx, "doi_normalised"] = fallback_doi
            if metrics is not None:
                metrics.record_application(
                    {
                        "pmid": pmid,
                        "fallback_doi": fallback_doi,
                        "previous_doi": existing_doi,
                        "sources": sources,
                        "action": "overwrite" if existing_doi else "fill",
                    }
                )
        return frame

    @classmethod
    def concat_frames_incrementally(
        cls,
        frames: Iterator[pd.DataFrame],
        *,
        batch_size: int,
    ) -> pd.DataFrame:
        """Concatenate ``frames`` in bounded batches to reduce peak memory usage."""

        combined: pd.DataFrame | None = None
        buffer: list[pd.DataFrame] = []

        def _flush(
            buffer_frames: list[pd.DataFrame], current: pd.DataFrame | None
        ) -> pd.DataFrame:
            """Concatenate ``buffer_frames`` into ``current`` and return the result."""

            if not buffer_frames:
                return (
                    current
                    if current is not None
                    else build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
                )

            if len(buffer_frames) == 1:
                batch = buffer_frames[0]
            else:
                batch = pd.concat(buffer_frames, ignore_index=True)

            if current is None:
                return batch

            return pd.concat([current, batch], ignore_index=True)

        for frame in frames:
            buffer.append(frame)
            if len(buffer) >= batch_size:
                combined = _flush(buffer, combined)
                buffer = []

        if combined is None and not buffer:
            return build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)

        return _flush(buffer, combined)

    def fetch_pubmed_records(
        self,
        pmids: Iterable[str],
        *,
        cfg: Config | None = None,
        sleep: float | None = None,
        pubmed_cfg: PubMedCfg | None = None,
        semantic_scholar_cfg: SemanticScholarCfg | None = None,
        openalex_cfg: OpenAlexCfg | None = None,
        crossref_cfg: CrossRefCfg | None = None,
        max_workers: int | None = None,
        batch_size: int | None = None,
        fallback_doi_map: Mapping[str, str] | None = None,
        return_generator: bool = False,
    ) -> pd.DataFrame | Iterator[pd.DataFrame]:
        """Retrieve metadata for a sequence of PubMed identifiers."""

        settings = cfg or self.cfg or Config()

        if sleep is None:
            sleep = settings.document.pubmed.sleep
        if pubmed_cfg is None:
            pubmed_cfg = settings.pubmed
        if semantic_scholar_cfg is None:
            semantic_scholar_cfg = settings.semantic_scholar
        if openalex_cfg is None:
            openalex_cfg = settings.openalex
        if crossref_cfg is None:
            crossref_cfg = settings.crossref
        if max_workers is None:
            max_workers = settings.document.pubmed.workers
        if batch_size is None:
            batch_size = settings.document.pubmed.batch_size

        openalex_limiter = get_limiter("openalex", openalex_cfg.rps, openalex_cfg.burst)
        crossref_limiter = get_limiter("crossref", crossref_cfg.rps, crossref_cfg.burst)

        rate_cfg = settings.rate
        system_limiter = None
        if (rate_cfg.global_rps or 0) > 0:
            system_limiter = get_global_limiter(
                rate_cfg.global_rps, rate_cfg.global_burst
            )

        def _service_limiter(
            name: str,
            *,
            rps: int | None,
            burst: int | None,
        ) -> RateLimiter | None:
            effective_rps = rps if rps is not None else rate_cfg.global_rps
            effective_burst = burst if burst is not None else rate_cfg.global_burst
            if effective_rps is None or effective_rps <= 0:
                return None
            if effective_burst is not None and effective_burst <= 0:
                effective_burst = None
            if (
                system_limiter is not None
                and (
                    (rps is None and burst is None)
                    or (
                        effective_rps == rate_cfg.global_rps
                        and effective_burst == rate_cfg.global_burst
                    )
                )
            ):
                return None
            return get_limiter(f"documents_{name}", effective_rps, effective_burst)

        pubmed_service_limiter = _service_limiter(
            "pubmed",
            rps=getattr(pubmed_cfg, "rps", None),
            burst=getattr(pubmed_cfg, "burst", None),
        )
        semantic_service_limiter = _service_limiter(
            "semantic_scholar",
            rps=getattr(semantic_scholar_cfg, "rps", None),
            burst=getattr(semantic_scholar_cfg, "burst", None),
        )
        openalex_service_limiter = _service_limiter(
            "openalex",
            rps=getattr(openalex_cfg, "rps", None),
            burst=getattr(openalex_cfg, "burst", None),
        )
        crossref_service_limiter = _service_limiter(
            "crossref",
            rps=getattr(crossref_cfg, "rps", None),
            burst=getattr(crossref_cfg, "burst", None),
        )

        def _acquire_documents(
            limiter: RateLimiter | None,
            *,
            use_global: bool = True,
            global_rate_limiter: RateLimiter | None = system_limiter,
        ) -> None:
            if use_global and global_rate_limiter is not None:
                global_rate_limiter.acquire()
            if limiter is not None:
                limiter.acquire()

        def _failure_records(
            batch: Sequence[str], message: str, *, source: str = "pubmed"
        ) -> list[dict[str, str]]:
            """Return placeholder rows describing a failure for ``batch`` PMIDs."""

            records: list[dict[str, str]] = []
            for pmid in batch:
                metadata = {"fetch_status": "error", "error_source": source}
                pubmed = {"PubMed.PMID": pmid, "PubMed.Error": message}
                scholar = {"scholar.PMID": pmid, "scholar.Error": message}
                openalex = {"OpenAlex.Error": message}
                crossref = {"crossref.Error": message}
                records.append(
                    merge_metadata(metadata, pubmed, scholar, openalex, crossref)
                )
            return records

        def _coerce_batch_argument(*candidates: object) -> list[str]:
            """Return the first iterable batch argument from ``candidates``."""

            for candidate in candidates:
                if isinstance(candidate, Sequence) and not isinstance(
                    candidate, (str, bytes, bytearray)
                ):
                    return [str(item) for item in candidate]
            raise TypeError(
                "fetch_pubmed_records() expected a sequence of PMIDs but received"
                f" {tuple(type(c).__name__ for c in candidates)}"
            )

        session_cfg = settings

        def _executor_capacity(limiter: RateLimiter | None, burst: int | None) -> int:
            limit = burst if burst is not None else 1
            if limiter is not None:
                limiter_burst = getattr(limiter, "burst", limit)
                try:
                    limiter_burst = int(limiter_burst)
                except (TypeError, ValueError):  # pragma: no cover - defensive
                    limiter_burst = limit
                limit = min(limit, limiter_burst)
            return max(1, limit)

        class _ThreadResources:
            def __init__(
                self,
                stack: ExitStack,
                base_session: requests.Session,
                session_pools: dict[str, _SessionPool],
            ) -> None:
                self.stack = stack
                self.base_session = base_session
                self.session_pools = session_pools
                self.sessions: dict[str, requests.Session] = {}
                self.session_finalizers: dict[str, Callable[[], None]] = {}

        class _SessionPool:
            def __init__(
                self,
                stack: ExitStack,
                factory: Callable[[], AbstractContextManager[requests.Session]],
                resources: _ThreadResources,
                service: str,
            ) -> None:
                self._stack = stack
                self._factory = factory
                self._resources = resources
                self._service = service

            @contextmanager
            def session(self) -> Iterator[requests.Session]:
                cached = self._resources.sessions.get(self._service)
                if cached is None:
                    context = self._factory()
                    session = context.__enter__()

                    def _finalizer(
                        _context: AbstractContextManager[requests.Session] = context,
                    ) -> None:
                        _context.__exit__(None, None, None)

                    self._resources.sessions[self._service] = session
                    self._resources.session_finalizers[self._service] = _finalizer
                    self._stack.callback(_finalizer)
                    cached = session
                yield cached

        thread_local_state = local()

        def _close_thread_resources(resources: _ThreadResources) -> None:
            try:
                resources.stack.close()
            except Exception:  # pragma: no cover - defensive
                logger.exception("thread_resources_close_failed")

        def _get_thread_resources() -> _ThreadResources:
            resources = getattr(thread_local_state, "resources", None)
            if resources is not None:
                return resources

            session_stack = ExitStack()
            base_session = session_stack.enter_context(
                session_with_retry(session_cfg.api, session_cfg.retry)
            )

            def _service_factory(
                mailto: str,
            ) -> Callable[[], AbstractContextManager[requests.Session]]:
                def _factory() -> AbstractContextManager[requests.Session]:
                    @contextmanager
                    def _context() -> Iterator[requests.Session]:
                        with session_with_retry(
                            session_cfg.api, session_cfg.retry
                        ) as derived_session:
                            if mailto and hasattr(derived_session, "headers"):
                                derived_session.headers["mailto"] = mailto
                            yield derived_session

                    return _context()

                return _factory

            session_factories: dict[
                str, Callable[[], AbstractContextManager[requests.Session]]
            ] = {
                "openalex": _service_factory(openalex_cfg.mailto),
                "crossref": _service_factory(crossref_cfg.mailto),
            }

            resources = _ThreadResources(session_stack, base_session, {})
            resources.session_pools = {
                service: _SessionPool(session_stack, factory, resources, service)
                for service, factory in session_factories.items()
            }
            setattr(thread_local_state, "resources", resources)

            finalizer = weakref.finalize(
                thread_local_state, _close_thread_resources, resources
            )
            setattr(thread_local_state, "resources_finalizer", finalizer)

            return resources

        def _fetch_batch(
            first: object,
            *rest: object,
            __cfg: Config = session_cfg,
            **__: object,
        ) -> list[dict[str, str]]:
            batch_list = _coerce_batch_argument(first, *rest)

            def _summarise_batch(pmids: Sequence[str]) -> dict[str, object]:
                sample_limit = 5
                sample = [pmids[index] for index in range(min(len(pmids), sample_limit))]
                if len(pmids) > sample_limit:
                    sample.append("...")
                return {"pmids_count": len(pmids), "pmids_sample": sample}

            batch_summary = _summarise_batch(batch_list)

            try:
                resources = _get_thread_resources()
                base_session = resources.base_session
                session_pools = resources.session_pools

                def _invoke_with_session(
                    factory: Callable[[], AbstractContextManager[requests.Session]],
                    fetcher: Callable[
                        [requests.Session, str, Any, RateLimiter | None], dict[str, str]
                    ],
                    identifier: str,
                    *,
                    cfg_obj: Any,
                    limiter: RateLimiter | None,
                ) -> dict[str, str]:
                    with factory() as nested_session:
                        return fetcher(nested_session, identifier, cfg_obj, limiter)

                _acquire_documents(pubmed_service_limiter)
                pubmed_list = pl.fetch_pubmed_batch(
                    base_session,
                    batch_list,
                    sleep,
                    cfg=pubmed_cfg,
                    retry_cfg=session_cfg.retry,
                )

                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]
                semantic_pmids = [pmid for pmid in pmids_in_batch if pmid]

                semsch_map: dict[str, dict[str, str]] = {}
                if semantic_pmids:
                    _acquire_documents(semantic_service_limiter)
                    semsch_list = ssl.fetch_semantic_scholar_batch(
                        base_session, semantic_pmids, sleep, cfg=semantic_scholar_cfg
                    )

                    semsch_map = {
                        s.get("scholar.PMID"): s
                        for s in semsch_list
                        if s.get("scholar.PMID")
                    }

                    fallback_pmids: list[str] = []
                    seen: set[str] = set()
                    for pmid in semantic_pmids:
                        if pmid in seen:
                            continue
                        seen.add(pmid)
                        record = semsch_map.get(pmid)
                        if record is None or record.get("scholar.Error"):
                            fallback_pmids.append(pmid)
                    for pmid in fallback_pmids:
                        _acquire_documents(semantic_service_limiter)
                        fallback_record = ssl.fetch_semantic_scholar(
                            base_session, pmid, sleep, cfg=semantic_scholar_cfg
                        )
                        semsch_map[pmid] = fallback_record

                combined_records: list[dict[str, str]] = []

                plan: list[tuple[int, dict[str, str], dict[str, str], str, str]] = []
                openalex_lookup: dict[str, list[int]] = {}
                crossref_lookup: dict[str, list[int]] = {}
                openalex_total = 0
                crossref_total = 0

                for index, pubmed in enumerate(pubmed_list):
                    pmid = pubmed.get("PubMed.PMID", "")
                    semsch = semsch_map.get(pmid, {}) if pmid else {}
                    fallback_doi = ""
                    if fallback_doi_map:
                        fallback_doi = fallback_doi_map.get(pmid, "")
                    doi = (
                        pubmed.get("PubMed.DOI")
                        or semsch.get("scholar.DOI")
                        or fallback_doi
                        or ""
                    )
                    plan.append((index, pubmed, semsch, pmid, doi))
                    if pmid:
                        openalex_lookup.setdefault(pmid, []).append(index)
                        openalex_total += 1
                    if doi:
                        crossref_lookup.setdefault(doi, []).append(index)
                        crossref_total += 1

                openalex_results: dict[int, dict[str, str]] = {}
                crossref_results: dict[int, dict[str, str]] = {}

                openalex_jobs = list(openalex_lookup.keys())

                def _fetch_openalex_job(pmid: str) -> dict[str, str]:
                    _acquire_documents(openalex_service_limiter, use_global=False)
                    return _invoke_with_session(
                        session_pools["openalex"].session,
                        ocl.fetch_openalex,
                        pmid,
                        cfg_obj=openalex_cfg,
                        limiter=openalex_limiter,
                    )

                if openalex_total:
                    logger.info(
                        "documents_cache_reuse",
                        service="openalex",
                        total=openalex_total,
                        unique=len(openalex_jobs),
                        hits=openalex_total - len(openalex_jobs),
                    )

                openalex_by_key: dict[str, dict[str, str]] = {}

                if openalex_jobs:
                    if openalex_executor is None:
                        for pmid in openalex_jobs:
                            openalex_by_key[pmid] = _fetch_openalex_job(pmid)
                    else:
                        future_to_key = {
                            openalex_executor.submit(_fetch_openalex_job, pmid): pmid
                            for pmid in openalex_jobs
                        }
                        for future in as_completed(future_to_key):
                            key = future_to_key[future]
                            openalex_by_key[key] = future.result()

                for pmid, indexes in openalex_lookup.items():
                    result = openalex_by_key.get(pmid, {})
                    for idx in indexes:
                        openalex_results[idx] = result

                crossref_jobs = [doi for doi in crossref_lookup.keys() if doi]

                def _fetch_crossref_job(doi: str) -> dict[str, str]:
                    _acquire_documents(crossref_service_limiter, use_global=False)
                    return _invoke_with_session(
                        session_pools["crossref"].session,
                        ocl.fetch_crossref,
                        doi,
                        cfg_obj=crossref_cfg,
                        limiter=crossref_limiter,
                    )

                if crossref_total:
                    logger.info(
                        "documents_cache_reuse",
                        service="crossref",
                        total=crossref_total,
                        unique=len(crossref_jobs),
                        hits=crossref_total - len(crossref_jobs),
                    )

                crossref_by_key: dict[str, dict[str, str]] = {}

                if crossref_jobs:
                    if crossref_executor is None:
                        for doi in crossref_jobs:
                            crossref_by_key[doi] = _fetch_crossref_job(doi)
                    else:
                        future_to_key = {
                            crossref_executor.submit(_fetch_crossref_job, doi): doi
                            for doi in crossref_jobs
                        }
                        for future in as_completed(future_to_key):
                            key = future_to_key[future]
                            crossref_by_key[key] = future.result()

                for doi, indexes in crossref_lookup.items():
                    result = crossref_by_key.get(doi, {})
                    for idx in indexes:
                        crossref_results[idx] = result

                for index, pubmed, semsch, pmid, _ in plan:
                    openalex = openalex_results.get(index, {}) if pmid else {}
                    crossref = crossref_results.get(index, {})
                    combined = merge_metadata(pubmed, semsch, openalex, crossref)
                    combined_records.append(combined)
                return combined_records
            except requests.RequestException as exc:  # pragma: no cover - network errors
                logger.warning(
                    "pubmed_batch_request_failed",
                    **batch_summary,
                    error=str(exc),
                )
                return _failure_records(batch_list, str(exc), source="pubmed")
            except Exception as exc:  # pragma: no cover - defensive safety net
                logger.warning(
                    "pubmed_batch_unexpected_error",
                    **batch_summary,
                    error=str(exc),
                )
                return _failure_records(batch_list, str(exc), source="pubmed")

        openalex_capacity = _executor_capacity(openalex_limiter, openalex_cfg.burst)
        crossref_capacity = _executor_capacity(crossref_limiter, crossref_cfg.burst)

        downstream_capacity = max(
            1,
            min(
                max_workers,
                openalex_capacity or max_workers,
                crossref_capacity or max_workers,
            ),
        )

        openalex_executor: ThreadPoolExecutor | None = None
        crossref_executor: ThreadPoolExecutor | None = None

        iterator = (p for p in pmids if p)

        tasks: dict[Future[list[dict[str, str]]], tuple[int, list[str]]] = {}
        frame_heap: list[tuple[int, pd.DataFrame]] = []
        next_to_emit = 0
        processed = 0
        completed_batches = 0
        max_in_flight = downstream_capacity * 2

        stack = ExitStack()
        batch_executor = stack.enter_context(
            ThreadPoolExecutor(max_workers=max_workers)
        )
        if openalex_capacity > 1:
            openalex_executor = stack.enter_context(
                ThreadPoolExecutor(max_workers=openalex_capacity)
            )
        if crossref_capacity > 1:
            crossref_executor = stack.enter_context(
                ThreadPoolExecutor(max_workers=crossref_capacity)
            )

        offset = 0
        pending: set[Future[list[dict[str, str]]]] = set()

        def _iter_pending(
            futures: Iterable[Future[list[dict[str, str]]]]
        ) -> Iterator[Future[list[dict[str, str]]]]:
            """Yield completed futures without copying the ``futures`` collection."""

            yield from as_completed(futures)

        lock = Lock()

        def _emit_ready_frames() -> Iterator[pd.DataFrame]:
            nonlocal next_to_emit

            with lock:
                while frame_heap and frame_heap[0][0] == next_to_emit:
                    _, frame = heapq.heappop(frame_heap)
                    yield frame
                    next_to_emit += len(frame)

        def _drain_future(
            done_future: Future[list[dict[str, str]]],
        ) -> Iterator[pd.DataFrame]:
            nonlocal processed, completed_batches

            pending.remove(done_future)
            batch_id, batch_pmids = tasks.pop(done_future)
            records = done_future.result()
            frame = build_dataframe(records, columns=DOCUMENT_SCHEMA_COLUMNS)
            with lock:
                heapq.heappush(frame_heap, (batch_id, frame))
            processed += len(batch_pmids)
            completed_batches += 1
            log_kwargs = {"count": processed, "batches": completed_batches}
            if completed_batches % self.DOCUMENT_PROGRESS_INFO_INTERVAL == 0:
                logger.info("documents_processed", **log_kwargs)
            else:
                logger.debug("documents_processed", **log_kwargs)

            yield from _emit_ready_frames()

        def _iter_frame_batches() -> Iterator[pd.DataFrame]:
            nonlocal offset

            for batch in _chunked(iterator, batch_size):
                while frame_heap and len(frame_heap) >= max_in_flight and pending:
                    done_future = next(as_completed(list(pending)))
                    yield from _drain_future(done_future)
                if not batch:
                    continue
                future = batch_executor.submit(_fetch_batch, batch)
                tasks[future] = (offset, batch)
                pending.add(future)
                offset += len(batch)
                if len(pending) >= max_in_flight:
                    done_future = next(_iter_pending(pending))
                    yield from _drain_future(done_future)

                yield from _emit_ready_frames()

            for done_future in _iter_pending(pending):
                yield from _drain_future(done_future)

            pending.clear()

            yield from _emit_ready_frames()

        def _iter_frames() -> Iterator[pd.DataFrame]:
            try:
                yield from _iter_frame_batches()
            finally:
                stack.close()
                finalizer = getattr(thread_local_state, "resources_finalizer", None)
                if finalizer is not None:
                    try:
                        finalizer()
                    finally:
                        delattr(thread_local_state, "resources_finalizer")
                if hasattr(thread_local_state, "resources"):
                    delattr(thread_local_state, "resources")

        frame_iter = _iter_frames()
        if return_generator:
            return frame_iter

        return self.concat_frames_incrementally(
            frame_iter,
            batch_size=self.DOCUMENT_FRAME_CONCAT_STRIDE,
        )


__all__ = [
    "DocumentPipeline",
    "FallbackDoiMetrics",
    "FallbackDoiState",
]
