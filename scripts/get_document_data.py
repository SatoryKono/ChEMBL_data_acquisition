"""Command line interface for retrieving document metadata from external sources.

The tool integrates :mod:`library.pubmed_library` and
:mod:`library.chembl_library` to collect information about publications from
several public APIs.  The interface mirrors :mod:`scripts.get_target_data` and
provides three sub-commands:

``pubmed``
    Query PubMed, Semantic Scholar, OpenAlex and CrossRef for a list of PMIDs.
``chembl``
    Retrieve document information from the ChEMBL API.
``all``
    Run the ChEMBL and PubMed pipelines and merge the results.

Example
-------
Fetch PubMed metadata for identifiers listed in ``pmids.csv``::

    python scripts/get_document_data.py pubmed --config config/config.yaml --input pmids.csv --output output.csv

The input file must contain a ``PMID`` column.

"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from contextlib import ExitStack, contextmanager, AbstractContextManager
from itertools import chain, islice, tee
from pathlib import Path
from threading import Lock
from typing import Any, cast

import pandas as pd
import requests
from pandera.errors import SchemaErrors


def _ensure_project_root() -> None:
    """Ensure the repository root is discoverable when executed as a script."""

    script_path = Path(__file__).resolve()
    project_root = script_path.parents[1]
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)


if __package__ in {None, ""}:
    _ensure_project_root()

from library import chembl_library as cl
from library import cli
from library import document_postprocessing as dp
from library import io
from library.csv_utils import write_csv_chunks_deterministic
from library import openalex_crossref_library as ocl
from library import pubmed_library as pl
from library import semantic_scholar_library as ssl
from library.clients import ChemblClient, _chunked
from library.cli import (
    LoggerConfig,
    build_root_parser,
    configure_logger,
    path_argument,
    positive_int,
    prepare_io_paths,
)
from library.config import (
    Config,
    CrossRefCfg,
    OpenAlexCfg,
    PubMedCfg,
    SemanticScholarCfg,
    _serialize_paths,
    ensure_dirs,
    openalex_session,
    print_config,
    session_with_retry,
    crossref_session,
)
from library.document_pipeline import (
    DOCUMENT_SCHEMA_COLUMNS,
    build_dataframe,
    build_quality_report,
    dataframe_to_strings,
    merge_metadata,
    merge_with_chembl,
    normalise_doi,
    save_quality_report,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.pipeline_metadata import add_pipeline_metadata
from library.rate_limiter import RateLimiter, get_limiter
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from schemas import DocumentsSchema, normalize_documents


DEFAULT_INPUT_NAME = "document.csv"
DEFAULT_OUTPUT_STEM = "documents"


def _build_fallback_doi_map(
    frame: pd.DataFrame,
    *,
    pmid_column: str,
    doi_column: str,
) -> dict[str, str]:
    """Return a mapping of PubMed identifiers to DOI overrides."""

    if pmid_column not in frame.columns or doi_column not in frame.columns:
        missing = [
            column
            for column in (pmid_column, doi_column)
            if column not in frame.columns
        ]
        raise ValueError(f"missing columns in fallback DOI file: {', '.join(missing)}")

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
    return mapping


def fetch_pubmed_records(
    pmids: Iterable[str],
    *args: object,
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
    """Retrieve metadata for a sequence of PubMed identifiers.

    Parameters
    ----------
    pmids : Iterable[str]
        Identifiers to query across the PubMed, Semantic Scholar, OpenAlex and
        CrossRef services.
    sleep : float, optional
        Seconds to pause between successive PubMed requests. The same interval
        is reused when the pipeline falls back to the Semantic Scholar
        single-record endpoint, preventing the thread pool from overwhelming
        either API.
    pubmed_cfg : PubMedCfg, optional
        Configuration used when requesting PubMed batches, including base URLs,
        authentication headers and retry behaviour.
    semantic_scholar_cfg : SemanticScholarCfg, optional
        Settings for the Semantic Scholar batch and fallback requests.
    openalex_cfg : OpenAlexCfg, optional
        Connection parameters for OpenAlex lookups.
    crossref_cfg : CrossRefCfg, optional
        Connection parameters for CrossRef lookups driven by the resolved DOI.
    max_workers : int, optional
        Maximum number of concurrent worker threads submitting batches. Higher
        values increase parallelism across both the PubMed batching and the
        Semantic Scholar enrichment stages.
    batch_size : int, optional
        Maximum number of PMIDs per PubMed request. Smaller batches reduce the
        risk of request failures at the cost of more round-trips and additional
        Semantic Scholar fallbacks.
    fallback_doi_map : Mapping[str, str], optional
        Explicit DOI overrides keyed by PMID. When provided, entries supply the
        DOI used to seed downstream CrossRef lookups whenever neither PubMed
        nor Semantic Scholar returns one.
    return_generator : bool, optional
        When ``True`` return a generator yielding ordered batches instead of a
        concatenated dataframe.

    Returns
    -------
    pandas.DataFrame or Iterator[pandas.DataFrame]
        Combined metadata from the different sources. When
        ``return_generator`` is ``True`` a generator yielding ordered
        :class:`~pandas.DataFrame` batches is returned instead of a single
        concatenated frame.

    Raises
    ------
    TypeError
        Raised when required configuration arguments are omitted or supplied
        multiple times. The function accepts either a
        :class:`~library.config.Config` instance or explicit keyword arguments
        describing the external services.

    Notes
    -----
    For backward compatibility the function also accepts a
    :class:`~library.config.Config` instance as the first positional argument
    after ``pmids``. When supplied, defaults are loaded from
    ``document.pubmed`` and related API sections before calling helper
    functions, and any keyword arguments supplied to
    :func:`fetch_pubmed_records` override those derived settings.

    """

    cfg: Config | None = None
    positional = list(args)
    if positional:
        candidate = positional[0]
        if isinstance(candidate, Config):
            cfg = candidate
            positional = positional[1:]
        elif sleep is None:
            if len(positional) < 4:
                raise TypeError(
                    "fetch_pubmed_records() missing required positional arguments: "
                    "'sleep', 'semantic_scholar_cfg', 'openalex_cfg', 'crossref_cfg'"
                )
            sleep = cast(float, positional[0])
            semantic_scholar_cfg = cast(SemanticScholarCfg, positional[1])
            openalex_cfg = cast(OpenAlexCfg, positional[2])
            crossref_cfg = cast(CrossRefCfg, positional[3])
            if len(positional) > 4:
                max_workers = cast(int, positional[4])
            if len(positional) > 5:
                batch_size = cast(int, positional[5])
            positional = []
        else:
            raise TypeError(
                "fetch_pubmed_records() received multiple values for 'sleep'"
            )
    if positional:
        raise TypeError("fetch_pubmed_records() got unexpected positional arguments")

    if cfg is not None:
        if sleep is None:
            sleep = cfg.document.pubmed.sleep
        if pubmed_cfg is None:
            pubmed_cfg = cfg.pubmed
        if semantic_scholar_cfg is None:
            semantic_scholar_cfg = cfg.semantic_scholar
        if openalex_cfg is None:
            openalex_cfg = cfg.openalex
        if crossref_cfg is None:
            crossref_cfg = cfg.crossref
        if max_workers is None:
            max_workers = cfg.document.pubmed.workers
        if batch_size is None:
            batch_size = cfg.document.pubmed.batch_size

    if (
        sleep is None
        or semantic_scholar_cfg is None
        or openalex_cfg is None
        or crossref_cfg is None
    ):
        raise TypeError(
            "fetch_pubmed_records() missing required configuration. "
            "Provide either a Config instance as the second positional argument "
            "or explicit keyword arguments."
        )

    if max_workers is None:
        max_workers = 1
    if batch_size is None:
        batch_size = 100

    openalex_limiter = get_limiter("openalex", openalex_cfg.rps, openalex_cfg.burst)
    crossref_limiter = get_limiter("crossref", crossref_cfg.rps, crossref_cfg.burst)

    settings = cfg or Config()
    rate_cfg = settings.rate
    if pubmed_cfg is None:
        pubmed_cfg = settings.pubmed

    documents_limiter = get_limiter(
        "documents_global", rate_cfg.global_rps, rate_cfg.global_burst
    )

    def _service_limiter(
        name: str,
        *,
        rps: int | None,
        burst: int | None,
    ) -> RateLimiter | None:
        effective_rps = rps if rps is not None else rate_cfg.global_rps
        effective_burst = burst if burst is not None else rate_cfg.global_burst
        if (
            rps is None
            and burst is None
            or (
                effective_rps == rate_cfg.global_rps
                and effective_burst == rate_cfg.global_burst
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
        limiter: RateLimiter | None, *, use_global: bool = True
    ) -> None:
        if use_global and documents_limiter is not None:
            documents_limiter.acquire()
        if limiter is not None:
            limiter.acquire()

    def _failure_records(batch: Sequence[str], message: str) -> list[dict[str, str]]:
        """Return placeholder rows describing a failure for ``batch`` PMIDs."""

        records: list[dict[str, str]] = []
        for pmid in batch:
            pubmed = {"PubMed.PMID": pmid, "PubMed.Error": message}
            scholar = {"scholar.PMID": pmid, "scholar.Error": message}
            openalex = {"OpenAlex.Error": message}
            crossref = {"crossref.Error": message}
            records.append(merge_metadata(pubmed, scholar, openalex, crossref))
        return records

    def _coerce_batch_argument(*candidates: object) -> list[str]:
        """Return the first iterable batch argument from ``candidates``."""

        for candidate in candidates:
            if isinstance(candidate, Sequence) and not isinstance(
                candidate, str | bytes | bytearray
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
            limit = min(limit, limiter.burst)
        return max(1, limit)

    def _fetch_batch(
        first: object,
        *rest: object,
        __cfg: Config = session_cfg,
        **__: object,
    ) -> list[dict[str, str]]:
        """Fetch metadata for a batch of PMIDs.

        Each worker opens its own :class:`requests.Session` and retrieves PubMed
        entries for all PMIDs in ``batch`` using a single request. Metadata from
        Semantic Scholar, OpenAlex and CrossRef are then fetched individually
        for each PMID. Exceptions are logged so a failure in one batch does not
        abort the whole process.
        """

        batch_list = _coerce_batch_argument(first, *rest)

        def _summarise_batch(pmids: Sequence[str]) -> dict[str, object]:
            sample_limit = 5
            sample = [pmids[index] for index in range(min(len(pmids), sample_limit))]
            if len(pmids) > sample_limit:
                sample.append("...")
            return {"pmids_count": len(pmids), "pmids_sample": sample}

        batch_summary = _summarise_batch(batch_list)

        class _SessionPool:
            def __init__(
                self,
                stack: ExitStack,
                factory: Callable[[], AbstractContextManager[requests.Session]],
            ) -> None:
                self._stack = stack
                self._factory = factory
                self._available: list[requests.Session] = []
                self._lock = Lock()

            def _acquire(self) -> requests.Session:
                with self._lock:
                    if self._available:
                        return self._available.pop()
                session_cm = self._factory()
                return self._stack.enter_context(session_cm)

            def _release(self, session: requests.Session) -> None:
                with self._lock:
                    self._available.append(session)

            @contextmanager
            def session(self) -> Iterator[requests.Session]:
                nested_session = self._acquire()
                try:
                    yield nested_session
                finally:
                    self._release(nested_session)

        try:
            with ExitStack() as session_stack:
                base_session = session_stack.enter_context(
                    session_with_retry(__cfg.api, __cfg.retry)
                )

                session_factories: dict[
                    str, Callable[[], AbstractContextManager[requests.Session]]
                ] = {
                    "openalex": lambda: openalex_session(
                        __cfg.api, __cfg.retry, openalex_cfg
                    ),
                    "crossref": lambda: crossref_session(
                        __cfg.api, __cfg.retry, crossref_cfg
                    ),
                }

                session_pools = {
                    service: _SessionPool(session_stack, factory)
                    for service, factory in session_factories.items()
                }

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
                    base_session, batch_list, sleep, cfg=pubmed_cfg
                )

                pmids_in_batch = [p.get("PubMed.PMID", "") for p in pubmed_list]
                semantic_pmids = [pmid for pmid in pmids_in_batch if pmid]

                semsch_map: dict[str, dict[str, str]] = {}
                if semantic_pmids:
                    # Fetch Semantic Scholar data in a single batch
                    _acquire_documents(semantic_service_limiter)
                    semsch_list = ssl.fetch_semantic_scholar_batch(
                        base_session, semantic_pmids, sleep, cfg=semantic_scholar_cfg
                    )

                    # Create a map for easy lookup
                    semsch_map = {
                        s.get("scholar.PMID"): s
                        for s in semsch_list
                        if s.get("scholar.PMID")
                    }

                    # Fallback to the single-record endpoint when the batch request fails
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
            return _failure_records(batch_list, str(exc))
        except Exception as exc:  # pragma: no cover - defensive safety net
            logger.warning(
                "pubmed_batch_unexpected_error",
                **batch_summary,
                error=str(exc),
            )
            return _failure_records(batch_list, str(exc))

    openalex_capacity = _executor_capacity(openalex_limiter, openalex_cfg.burst)
    crossref_capacity = _executor_capacity(crossref_limiter, crossref_cfg.burst)

    openalex_executor: ThreadPoolExecutor | None = None
    crossref_executor: ThreadPoolExecutor | None = None

    iterator = (p for p in pmids if p)

    tasks: dict[Future[list[dict[str, str]]], tuple[int, list[str]]] = {}
    completed: dict[int, list[dict[str, str]]] = {}
    next_to_emit = 0
    processed = 0
    max_in_flight = max(1, max_workers * 2)

    stack = ExitStack()
    batch_executor = stack.enter_context(ThreadPoolExecutor(max_workers=max_workers))
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

    def _drain_future(
        done_future: Future[list[dict[str, str]]],
    ) -> Iterator[list[dict[str, str]]]:
        nonlocal processed, next_to_emit

        pending.remove(done_future)
        batch_id, batch_pmids = tasks.pop(done_future)
        completed[batch_id] = done_future.result()
        processed += len(batch_pmids)
        logger.info("documents_processed", count=processed)

        while next_to_emit in completed:
            records = completed.pop(next_to_emit)
            yield records
            next_to_emit += len(records)

    def _emit_ready_batches() -> Iterator[list[dict[str, str]]]:
        nonlocal next_to_emit

        while next_to_emit in completed:
            records = completed.pop(next_to_emit)
            yield records
            next_to_emit += len(records)

    def _iter_records() -> Iterator[list[dict[str, str]]]:
        nonlocal offset

        for batch in _chunked(iterator, batch_size):
            if not batch:
                continue
            future = batch_executor.submit(_fetch_batch, batch)
            tasks[future] = (offset, batch)
            pending.add(future)
            offset += len(batch)
            if len(pending) >= max_in_flight:

                done_future = next(as_completed(list(pending)))
                yield from _drain_future(done_future)

            yield from _emit_ready_batches()

        for done_future in as_completed(list(pending)):
            yield from _drain_future(done_future)

        pending.clear()

        yield from _emit_ready_batches()

    def _iter_frames() -> Iterator[pd.DataFrame]:
        try:
            for records_batch in _iter_records():
                if not records_batch:
                    yield build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
                    continue
                yield build_dataframe(records_batch, columns=DOCUMENT_SCHEMA_COLUMNS)
        finally:
            stack.close()

    frame_iter = _iter_frames()
    if return_generator:
        return frame_iter

    frame_iter, concat_iter = tee(frame_iter)
    try:
        first_frame = next(concat_iter)
    except StopIteration:
        return build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
    return pd.concat(chain([first_frame], concat_iter), ignore_index=True)


_NUMERIC_EXPORT_COLUMNS = {
    "publication_type_score_review",
    "publication_type_score_experimental",
    "publication_type_score_unknown",
}


_EXPORT_COLUMNS = [
    "PubMed.PMID",
    "PubMed.DOI",
    "PubMed.ArticleTitle",
    "PubMed.Abstract",
    "PubMed.JournalTitle",
    "PubMed.JournalISOAbbrev",
    "PubMed.Volume",
    "PubMed.Issue",
    "PubMed.StartPage",
    "PubMed.EndPage",
    "PubMed.ISSN",
    "PubMed.PublicationType",
    "PubMed.MeSH_Descriptors",
    "PubMed.MeSH_Qualifiers",
    "PubMed.ChemicalList",
    "PubMed.YearCompleted",
    "PubMed.MonthCompleted",
    "PubMed.DayCompleted",
    "PubMed.YearRevised",
    "PubMed.MonthRevised",
    "PubMed.DayRevised",
    "PubMed.Error",
    "scholar.PMID",
    "scholar.DOI",
    "scholar.PublicationTypes",
    "scholar.Venue",
    "scholar.SemanticScholarId",
    "scholar.ExternalIds",
    "scholar.Error",
    "OpenAlex.PMID",
    "OpenAlex.DOI",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
    "OpenAlex.Genre",
    "OpenAlex.Venue",
    "OpenAlex.MeshDescriptors",
    "OpenAlex.MeshQualifiers",
    "OpenAlex.Id",
    "OpenAlex.Error",
    "crossref.DOI",
    "crossref.Type",
    "crossref.Subtype",
    "crossref.Title",
    "crossref.Subtitle",
    "crossref.Subject",
    "crossref.Error",
    "publication_types_normalised",
    "publication_review_score",
    "publication_experimental_score",
    "publication_class",
    "ChEMBL.document_chembl_id",
    "ChEMBL.title",
    "ChEMBL.abstract",
    "ChEMBL.doi",
    "ChEMBL.year",
    "ChEMBL.journal",
    "ChEMBL.journal_abbrev",
    "ChEMBL.volume",
    "ChEMBL.issue",
    "ChEMBL.first_page",
    "ChEMBL.last_page",
    "ChEMBL.pubmed_id",
    "ChEMBL.authors",
    "ChEMBL.source",
]

_EXPORT_COLUMN_RENAMES = {
    "document_chembl_id": "ChEMBL.document_chembl_id",
    "title": "ChEMBL.title",
    "abstract": "ChEMBL.abstract",
    "doi": "ChEMBL.doi",
    "year": "ChEMBL.year",
    "journal": "ChEMBL.journal",
    "journal_abbrev": "ChEMBL.journal_abbrev",
    "volume": "ChEMBL.volume",
    "issue": "ChEMBL.issue",
    "first_page": "ChEMBL.first_page",
    "last_page": "ChEMBL.last_page",
    "pubmed_id": "ChEMBL.pubmed_id",
    "authors": "ChEMBL.authors",
    "source": "ChEMBL.source",
    "publication_type_score_review": "publication_review_score",
    "publication_type_score_experimental": "publication_experimental_score",
}

_EXPORT_COALESCE_SOURCES = {
    "OpenAlex.PMID": ["OpenAlex.PMID", "PubMed.PMID", "scholar.PMID"],
    "OpenAlex.DOI": ["OpenAlex.DOI", "PubMed.DOI", "scholar.DOI", "doi_normalised"],
    "crossref.DOI": ["crossref.DOI", "doi_normalised", "PubMed.DOI", "scholar.DOI"],
}

_EXPORT_SORT_FALLBACK = [
    "ChEMBL.document_chembl_id",
    "PubMed.PMID",
    "scholar.PMID",
    "OpenAlex.PMID",
    "ChEMBL.pubmed_id",
]


_EXPORT_STREAM_CHUNK_SIZE = 10_000


def _iter_export_chunks(df: pd.DataFrame, *, chunk_size: int) -> Iterable[pd.DataFrame]:
    """Yield export-ready DataFrame chunks from ``df``."""

    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    if df.empty:
        yield build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False)
        return

    total = len(df)
    for start in range(0, total, chunk_size):
        stop = start + chunk_size
        chunk = df.iloc[start:stop]
        export_chunk = build_dataframe(
            chunk, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
        )
        export_chunk = dataframe_to_strings(export_chunk, skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(export_chunk)


def _coalesce_columns(df: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    """Return the first non-empty value across ``columns`` for each row."""

    result = pd.Series("", index=df.index, dtype=object)
    for col in columns:
        if col not in df.columns:
            continue
        values = df[col].fillna("").astype(str)
        mask = result.eq("")
        if mask.any():
            result.loc[mask] = values.loc[mask]
    return result


def _prepare_export_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Rename and project columns to match the export schema."""

    frame = df.copy()
    rename_map = {k: v for k, v in _EXPORT_COLUMN_RENAMES.items() if k in frame.columns}
    if rename_map:
        frame = frame.rename(columns=rename_map)

    for target, sources in _EXPORT_COALESCE_SOURCES.items():
        frame[target] = _coalesce_columns(frame, sources)

    for column in _EXPORT_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""

    return frame[_EXPORT_COLUMNS]


def _iter_export_chunks(
    df: pd.DataFrame,
    *,
    chunk_size: int | None,
) -> Iterable[pd.DataFrame]:
    """Yield export-ready chunks with deterministic column ordering."""

    total_rows = len(df)
    if total_rows == 0:
        empty = dataframe_to_strings(df.copy(), skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(empty)
        return

    effective_size = chunk_size if chunk_size and chunk_size > 0 else total_rows
    for start in range(0, total_rows, effective_size):
        stop = start + effective_size
        chunk = df.iloc[start:stop].copy()
        chunk = dataframe_to_strings(chunk, skip=_NUMERIC_EXPORT_COLUMNS)
        yield _prepare_export_frame(chunk)


def _resolve_chunk_size(value: int | None) -> int | None:
    """Return ``value`` when positive, otherwise ``None``."""

    if value is None:
        return None
    if value <= 0:
        logger.warning("invalid_csv_chunksize", value=value)
        return None
    return value


def _write_export_chunks(
    chunks: Iterable[pd.DataFrame],
    path: Path,
    *,
    cfg: Config,
    key_cols: Sequence[str],
    chunk_size: int | None,
) -> Path:
    """Stream ``chunks`` to ``path`` using the deterministic CSV writer."""

    if chunk_size:
        return write_csv_chunks_deterministic(
            chunks,
            path,
            col_order=list(_EXPORT_COLUMNS),
            key_cols=list(key_cols),
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            cfg=cfg,
            chunksize=chunk_size,
            merge_chunksize=chunk_size,
        )

    return write_csv_chunks_deterministic(
        chunks,
        path,
        col_order=list(_EXPORT_COLUMNS),
        key_cols=list(key_cols),
        sep=cfg.io.csv_sep,
        encoding=cfg.io.csv_encoding,
        cfg=cfg,
    )


def _load_export_ready_frame(path: Path, *, cfg: Config) -> pd.DataFrame:
    """Load the streamed CSV export for quality reporting."""

    try:
        loaded = pd.read_csv(
            path,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    except pd.errors.EmptyDataError:
        loaded = pd.DataFrame(columns=list(_EXPORT_COLUMNS))
    missing = [column for column in _EXPORT_COLUMNS if column not in loaded.columns]
    for column in missing:
        loaded[column] = ""
    return loaded[_EXPORT_COLUMNS]


def _finalise_export(
    df: pd.DataFrame | Iterable[pd.DataFrame],
    output: Path,
    cfg: Config,
    *,
    input_csv: Path,
    key_columns: Sequence[str] | None = None,
    chunk_size: int | None = None,
) -> int:
    """Validate input frames and write CSV/metadata artefacts."""

    if isinstance(df, pd.DataFrame):
        frames_iterable: Iterable[pd.DataFrame] = (df,)
    else:
        frames_iterable = df

    frames_iterator = iter(frames_iterable)
    required_cols = {
        name for name, col in DocumentsSchema.columns.items() if col.required
    }
    optional_cols = set(DocumentsSchema.columns) - required_cols

    present_columns: set[str] = set()
    missing_required: set[str] = set(required_cols)
    missing_optional: set[str] = set(optional_cols)

    stream_chunk = max(1, int(chunk_size or _EXPORT_STREAM_CHUNK_SIZE))
    failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
    errors = SidecarErrors()
    rows_total = 0
    rows_kept = 0
    exit_code = 0
    emitted_chunk = False

    def _validated_chunks() -> Iterator[pd.DataFrame]:
        nonlocal rows_total, rows_kept, exit_code, emitted_chunk
        nonlocal missing_required, missing_optional

        def _update_column_sets(df: pd.DataFrame) -> None:
            nonlocal missing_required, missing_optional
            present_columns.update(df.columns)
            missing_required = required_cols - present_columns
            missing_optional = optional_cols - present_columns

        for frame in frames_iterator:
            with_metadata = add_pipeline_metadata(frame)
            ordered = build_dataframe(
                with_metadata, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            _update_column_sets(ordered)
            rows_total += len(ordered)
            validated = ordered
            try:
                validated = DocumentsSchema.validate(ordered, lazy=True)
            except SchemaErrors as exc:
                for row in exc.failure_cases.to_dict("records"):
                    errors.add_error(row)
                logger.error(
                    "document_validation_failed",
                    failure_count=len(exc.failure_cases),
                    failure_path=str(failure_path),
                    error=str(exc),
                )
                validated = getattr(exc, "validated_data", ordered)
                exit_code = 1
            rows_kept += len(validated)
            cleaned = build_dataframe(
                validated, columns=DOCUMENT_SCHEMA_COLUMNS, fill_missing=False
            )
            for chunk in _iter_export_chunks(cleaned, chunk_size=stream_chunk):
                emitted_chunk = True
                yield chunk

        if not emitted_chunk:
            empty = build_dataframe(
                add_pipeline_metadata(pd.DataFrame()),
                columns=DOCUMENT_SCHEMA_COLUMNS,
                fill_missing=False,
            )
            _update_column_sets(empty)
            for chunk in _iter_export_chunks(empty, chunk_size=stream_chunk):
                emitted_chunk = True
                yield chunk

    export_generator = _validated_chunks()

    key_cols: list[str] = []
    if key_columns:
        for column in key_columns:
            mapped = _EXPORT_COLUMN_RENAMES.get(column, column)
            if mapped in _EXPORT_COLUMNS and mapped not in key_cols:
                key_cols.append(mapped)
    if not key_cols:
        for candidate in _EXPORT_SORT_FALLBACK:
            if candidate in _EXPORT_COLUMNS:
                key_cols = [candidate]
                break
    if not key_cols:
        key_cols = [_EXPORT_COLUMNS[0]]

    col_order = list(_EXPORT_COLUMNS)
    try:
        csv_path = write_csv_chunks_deterministic(
            export_generator,
            output,
            cfg=cfg,
            key_cols=key_cols,
            col_order=col_order,
            chunksize=stream_chunk,
            merge_chunksize=stream_chunk,
            sort_chunksize=stream_chunk,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
    except OSError as exc:
        logger.error("csv_write_failed", error=str(exc), path=str(output))
        return 1

    if missing_required:
        logger.warning(
            "validation_skipped_missing_required",
            columns=sorted(missing_required),
        )
        exit_code = 1
    elif missing_optional:
        logger.warning(
            "missing_optional_columns",
            columns=sorted(missing_optional),
        )

    errors.save(failure_path)

    rows_dropped = rows_total - rows_kept
    logger.info("write_done", rows=rows_kept, path=str(csv_path))

    try:
        export_ready = _load_export_ready_frame(csv_path, cfg=cfg)
    except (OSError, ValueError) as exc:
        logger.error("quality_frame_load_failed", error=str(exc), path=str(csv_path))
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(input_csv)},
        stats=stats,
        schema="DocumentsSchema",
    )

    quality_path = csv_path.with_suffix(".quality.json")
    try:
        report = build_quality_report(export_ready)
        save_quality_report(report, quality_path)
    except (OSError, TypeError, ValueError) as exc:
        logger.error(
            "quality_report_write_failed",
            error=str(exc),
            path=str(quality_path),
        )
        return 1

    try:
        analyze_table_quality(export_ready, table_name=str(csv_path.with_suffix("")))
    except ValueError as exc:
        logger.error("quality_report_generation_failed", error=str(exc))
        return 1
    return exit_code


def run_pubmed(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``pubmed`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration containing rate limiting, API and CSV export
        settings.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. A non-zero value indicates that an error occurred
        while reading input identifiers, fetching metadata or writing the
        resulting CSV.
    """
    pubmed_defaults = cfg.document.pubmed
    limit = getattr(args, "limit", pubmed_defaults.limit)
    if limit is not None and limit < 0:
        logger.error(
            "invalid_limit",
            section="document.pubmed",
            limit=limit,
        )
        return 1
    try:
        pmids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", pubmed_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1
    offset = getattr(args, "offset", 0)
    if offset:
        pmids_iter = islice(pmids_iter, offset, None)
        logger.info("process_offset", offset=offset)
    pmids: Iterable[str] = pmids_iter
    if limit is not None:
        limited_pmids = list(islice(pmids_iter, limit))
        pmids = limited_pmids
        logger.info("process_limit", limit=len(limited_pmids))

    fallback_doi_map: Mapping[str, str] | None = None
    fallback_csv = getattr(args, "fallback_doi_csv", None)
    if fallback_csv:
        try:
            fallback_frame = io.read_csv(fallback_csv, cfg=cfg.io)
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "fallback_doi_read_failed",
                error=str(exc),
                path=str(fallback_csv),
            )
            return 1
        try:
            fallback_map = _build_fallback_doi_map(
                fallback_frame,
                pmid_column=getattr(args, "fallback_doi_pmid_column", "PMID"),
                doi_column=getattr(args, "fallback_doi_value_column", "DOI"),
            )
        except ValueError as exc:
            logger.error(
                "fallback_doi_invalid",
                error=str(exc),
                path=str(fallback_csv),
            )
            return 1
        fallback_doi_map = fallback_map or None

    try:
        frame_iter = fetch_pubmed_records(
            pmids,
            cfg,
            sleep=getattr(args, "sleep", pubmed_defaults.sleep),
            pubmed_cfg=cfg.pubmed,
            semantic_scholar_cfg=cfg.semantic_scholar,
            openalex_cfg=cfg.openalex,
            crossref_cfg=cfg.crossref,
            max_workers=getattr(args, "workers", pubmed_defaults.workers),
            batch_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
            fallback_doi_map=fallback_doi_map,
            return_generator=True,
        )
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
        normalised_frames = (normalize_documents(frame) for frame in frame_iter)
        exit_code = _finalise_export(
            normalised_frames,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "batch_size", pubmed_defaults.batch_size),
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("pubmed_pipeline_failed", error=str(exc))
        return 1
    return exit_code


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration providing ChEMBL client, retry and CSV export
        options.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. A non-zero value indicates that reading the input
        identifiers, fetching ChEMBL data or exporting the CSV failed. Network
        errors are logged and converted into placeholder rows where possible.
    """
    chembl_defaults = cfg.document.chembl
    limit = getattr(args, "limit", chembl_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.chembl", limit=limit)
        return 1

    # Configure session for ChEMBL requests
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids_iter = io.read_ids(
                args.input_csv,
                column=getattr(args, "column", chembl_defaults.column),
                cfg=cfg.io,
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error(
                "input_read_failed",
                error=str(exc),
                path=str(args.input_csv),
            )
            return 1

        offset = getattr(args, "offset", 0)
        if offset:
            ids_iter = islice(ids_iter, offset, None)
            logger.info("process_offset", offset=offset)

        ids: Iterable[str] = ids_iter
        if limit is not None:
            limited_ids = list(islice(ids_iter, limit))
            ids = limited_ids
            logger.info("process_limit", limit=len(limited_ids))

        try:
            df = cl.get_documents(
                ids,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
                timeout=getattr(args, "timeout", chembl_defaults.timeout),
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error(
                "chembl_documents_fetch_failed",
                error=str(exc),
                chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
            )
            return 1
        if "doi" in df.columns:
            df["doi"] = df["doi"].map(normalise_doi)
        output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
        df = normalize_documents(df)
        exit_code = _finalise_export(
            df,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", chembl_defaults.chunk_size),
        )
        return exit_code


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run ChEMBL and PubMed pipelines and merge their outputs.

    Parameters
    ----------
    cfg : Config
        Application configuration combining all document pipeline defaults.
    args : argparse.Namespace
        Parsed command-line arguments produced by :func:`build_parser`.

    Returns
    -------
    int
        ``0`` on success. Non-zero results indicate that reading identifiers,
        fetching data from upstream APIs or writing derived artefacts failed.

    Raises
    ------
    ValueError
        Raised when DOI fallback information derived from the ChEMBL payload is
        internally inconsistent.
    """
    all_defaults = cfg.document.all
    limit = getattr(args, "limit", all_defaults.limit)
    if limit is not None and limit < 0:
        logger.error("invalid_limit", section="document.all", limit=limit)
        return 1

    # Prepare shared session before performing any API calls
    try:
        ids_iter = io.read_ids(
            args.input_csv,
            column=getattr(args, "column", all_defaults.column),
            cfg=cfg.io,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(
            "input_read_failed",
            error=str(exc),
            path=str(args.input_csv),
        )
        return 1

    offset = getattr(args, "offset", 0)
    if offset:
        ids_iter = islice(ids_iter, offset, None)
        logger.info("process_offset", offset=offset)

    ids_source: Iterable[str] = ids_iter
    if limit is not None:
        limited_ids = list(islice(ids_iter, limit))
        ids_source = limited_ids
        logger.info("process_limit", limit=len(limited_ids))

    iterator = iter(ids_source)
    sample_size = getattr(args, "chunk_size", all_defaults.chunk_size)
    sample_ids = list(islice(iterator, sample_size))
    ids_for_fetch = chain(sample_ids, iterator)
    try:
        with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
            doc_df = cl.get_documents(
                ids_for_fetch,
                cfg=cfg.api,
                client=client,
                chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
                timeout=getattr(args, "timeout", all_defaults.timeout),
            )
    except (requests.RequestException, ValueError) as exc:
        logger.error(
            "chembl_documents_fetch_failed",
            ids=sample_ids,
            error=str(exc),
            chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
        )
        return 1
    output = Path(args.output_csv or io.default_output_path(args.input_csv, cfg.io))
    if "doi" in doc_df.columns:
        doc_df["doi"] = doc_df["doi"].map(normalise_doi)
    if doc_df.empty or "pubmed_id" not in doc_df:
        processed = dp.postprocess_documents(doc_df)
        extra_cols = [c for c in doc_df.columns if c not in processed.columns]
        if extra_cols:
            processed = processed.merge(
                doc_df[["document_chembl_id"] + extra_cols],
                on="document_chembl_id",
                how="left",
            )
        processed = normalize_documents(processed)
        exit_code = _finalise_export(
            processed,
            output,
            cfg,
            input_csv=Path(args.input_csv),
            key_columns=["document_chembl_id"],
            chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
        )
        return exit_code

    # Normalise PubMed identifiers to strings to avoid dtype mismatches
    pubmed_ids = pd.to_numeric(doc_df["pubmed_id"], errors="coerce").astype("Int64")
    doi_fallback_map: dict[str, str] = {}
    if "doi" in doc_df.columns:
        doi_series = doc_df["doi"].astype("string")
        mask = pubmed_ids.notna() & doi_series.notna() & (doi_series != "")
        if mask.any():
            masked_pmids = pubmed_ids[mask].tolist()
            masked_dois = doi_series[mask].tolist()
            if len(masked_pmids) != len(masked_dois):
                raise ValueError("mismatched DOI fallback map lengths")
            doi_fallback_map = {
                str(pmid): str(doi)
                for pmid, doi in zip(masked_pmids, masked_dois, strict=True)
            }
    pmids = pubmed_ids.dropna().astype(str).tolist()
    pubmed_frames = fetch_pubmed_records(
        pmids,
        cfg,
        sleep=getattr(args, "sleep", all_defaults.sleep),
        semantic_scholar_cfg=cfg.semantic_scholar,
        openalex_cfg=cfg.openalex,
        crossref_cfg=cfg.crossref,
        max_workers=getattr(args, "workers", all_defaults.workers),
        batch_size=getattr(args, "batch_size", all_defaults.batch_size),
        pubmed_cfg=cfg.pubmed,
        fallback_doi_map=doi_fallback_map or None,
        return_generator=True,
    )
    concat_iter = iter(pubmed_frames)
    try:
        first_frame = next(concat_iter)
    except StopIteration:
        pub_df = build_dataframe([], columns=DOCUMENT_SCHEMA_COLUMNS)
    else:
        pub_df = pd.concat(chain([first_frame], concat_iter), ignore_index=True)
    doc_df["pubmed_id"] = pubmed_ids.astype("Int64").astype("string").fillna("")
    merged = merge_with_chembl(doc_df, pub_df)
    processed = dp.postprocess_documents(merged)
    extra_cols = [c for c in merged.columns if c not in processed.columns]
    if extra_cols:
        processed = processed.merge(
            merged[["document_chembl_id"] + extra_cols],
            on="document_chembl_id",
            how="left",
        )
    processed = normalize_documents(processed)
    exit_code = _finalise_export(
        processed,
        output,
        cfg,
        input_csv=Path(args.input_csv),
        key_columns=["document_chembl_id"],
        chunk_size=getattr(args, "chunk_size", all_defaults.chunk_size),
    )
    return exit_code


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the argument parser for document utilities.

    Returns
    -------
    tuple[argparse.ArgumentParser, LoggerConfig]
        Parser populated with all sub-commands and default logging
        configuration used by :func:`main`.
    """
    root, shared, log_cfg = build_root_parser()
    root.set_defaults(input_csv=Path(DEFAULT_INPUT_NAME))
    parser = argparse.ArgumentParser(
        description="Document data utilities", parents=[root]
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pubmed = sub.add_parser(
        "pubmed", parents=[shared], help="Fetch data from PubMed and related APIs"
    )
    pubmed.add_argument(
        "--column", default="PMID", help="Column name containing identifiers"
    )
    pubmed.add_argument(
        "--sleep", type=float, default=5.0, help="Seconds to sleep between requests"
    )
    pubmed.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent requests"
    )
    pubmed.add_argument(
        "--batch-size",
        type=positive_int,
        default=100,
        help="Maximum PMIDs per PubMed request",
    )
    pubmed.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    pubmed.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    pubmed.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex",
    )
    pubmed.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef",
    )
    pubmed.add_argument(
        "--fallback-doi-csv",
        type=path_argument,
        default=None,
        help="Optional CSV file providing PMID to DOI overrides",
    )
    pubmed.add_argument(
        "--fallback-doi-pmid-column",
        default="PMID",
        help="Column containing PubMed identifiers in fallback CSV",
    )
    pubmed.add_argument(
        "--fallback-doi-value-column",
        default="DOI",
        help="Column containing DOI values in fallback CSV",
    )
    pubmed.set_defaults(func=run_pubmed)

    chembl = sub.add_parser(
        "chembl", parents=[shared], help="Fetch document information from ChEMBL"
    )
    chembl.add_argument(
        "--column",
        default="document_chembl_id",
        help="Column name containing identifiers",
    )
    chembl.add_argument(
        "--chunk-size",
        type=positive_int,
        default=5,
        help="Maximum number of IDs per request",
    )
    chembl.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    chembl.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    chembl.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    chembl.set_defaults(func=run_chembl)

    all_cmd = sub.add_parser(
        "all", parents=[shared], help="Run both ChEMBL and PubMed pipelines"
    )
    all_cmd.add_argument(
        "--column",
        default="document_chembl_id",
        help="Column in the input CSV",
    )
    all_cmd.add_argument(
        "--chunk-size",
        type=positive_int,
        default=5,
        help="Maximum IDs per request",
    )
    all_cmd.add_argument(
        "--sleep",
        type=float,
        default=5.0,
        help="Seconds to sleep between PubMed requests",
    )
    all_cmd.add_argument(
        "--workers", type=int, default=1, help="Number of concurrent PubMed requests"
    )
    all_cmd.add_argument(
        "--batch-size",
        type=positive_int,
        default=50,
        help="Maximum PMIDs per PubMed request",
    )
    all_cmd.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    all_cmd.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of identifiers to process",
    )
    all_cmd.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Number of identifiers to skip before processing",
    )
    all_cmd.add_argument(
        "--openalex-rps",
        type=float,
        default=None,
        help="Requests per second limit for OpenAlex",
    )
    all_cmd.add_argument(
        "--crossref-rps",
        type=float,
        default=None,
        help="Requests per second limit for CrossRef",
    )
    all_cmd.set_defaults(func=run_all)

    parser.subparsers_map = {  # type: ignore[attr-defined]
        "pubmed": pubmed,
        "chembl": chembl,
        "all": all_cmd,
    }

    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Parameters
    ----------
    argv : Sequence[str] | None, optional
        Command-line arguments to parse. When ``None`` the process arguments
        from :data:`sys.argv` are used.

    Returns
    -------
    int
        ``0`` when the selected pipeline completes successfully, non-zero
        otherwise.

    Raises
    ------
    SystemExit
        Raised when argument validation fails and ``argparse`` aborts
        execution.
    """
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    prepare_io_paths(
        args,
        input_default=DEFAULT_INPUT_NAME,
        output_stem=DEFAULT_OUTPUT_STEM,
    )
    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    limit_value = getattr(args, "limit", None)
    if limit_value is not None and limit_value <= 0:
        subparser.error("--limit must be a positive integer")
    offset_value = getattr(args, "offset", 0)
    if offset_value < 0:
        subparser.error("--offset must be zero or a positive integer")
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    mapping = {
        "column": f"document.{args.command}.column",
        "limit": f"document.{args.command}.limit",
    }
    if args.command == "pubmed":
        mapping.update(
            {
                "sleep": "document.pubmed.sleep",
                "workers": "document.pubmed.workers",
                "batch_size": "document.pubmed.batch_size",
            }
        )
    elif args.command == "chembl":
        mapping.update(
            {
                "chunk_size": "document.chembl.chunk_size",
                "timeout": "document.chembl.timeout",
            }
        )
    elif args.command == "all":
        mapping.update(
            {
                "chunk_size": "document.all.chunk_size",
                "sleep": "document.all.sleep",
                "workers": "document.all.workers",
                "batch_size": "document.all.batch_size",
                "timeout": "document.all.timeout",
            }
        )
    try:
        cfg: Config = cli.apply_config_overrides(
            args,
            subparser,
            args.config,
            mapping=mapping
            | {
                "openalex_rps": "openalex.rps",
                "crossref_rps": "crossref.rps",
            },
            base_parser=parser,
        )
        cfg.api.timeout_read = getattr(args, "timeout", cfg.api.timeout_read)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        output_path = Path(
            args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        )
        args.output_csv = output_path
        if args.skip_existing and output_path.exists() and not args.force:
            logger.info("pipeline_skip_existing", output=str(output_path))
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error("config_override_error", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("directory_setup_failed", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code = cast(int, args.func(cfg, args))
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
