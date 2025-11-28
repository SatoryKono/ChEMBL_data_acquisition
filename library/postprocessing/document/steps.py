"""Document acquisition helper steps combining ChEMBL and DOI metadata."""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass
from typing import Callable, Iterable, Mapping, MutableMapping, Sequence
from urllib.parse import parse_qsl, urlencode, urljoin, urlsplit, urlunsplit

import pandas as pd
from requests import Session

from library.clients import ChemblClient
from library.crossref.clients import crossref as crossref_client
from library.clients import openalex as openalex_client
from library.common.rate_limiter import RateLimiter, get_limiter
from library.config import Config, crossref_session, openalex_session
from library.schemas.document_schema import validate_document_frame
from library.utils.qc_report import TableQualityProfiler, build_reports_from_profiler

_CHEMBL_DOCUMENT_COLUMNS: tuple[str, ...] = (
    "document_chembl_id",
    "title",
    "doi",
    "pubmed_id",
    "doc_type",
    "journal",
    "year",
)

_CROSSREF_COLUMNS: tuple[str, ...] = (
    "doi_key",
    "crossref_title",
    "crossref_doc_type",
    "crossref_subject",
    "crossref_error",
)

_OPENALEX_COLUMNS: tuple[str, ...] = (
    "doi_key",
    "openalex_title",
    "openalex_doc_type",
    "openalex_type_crossref",
    "openalex_publication_year",
    "openalex_error",
)

_CHEMBL_PAGE_SIZE_MAX = 1000

_DOI_PREFIX_PATTERN = re.compile(r"^doi:\s*", flags=re.IGNORECASE)
_HTTP_DOI_PREFIX_PATTERN = re.compile(r"^https?://(?:dx\.)?doi\.org/")


@dataclass(frozen=True)
class DocumentFetchResult:
    """Container holding the enriched document frame and QC artefacts."""

    data: pd.DataFrame
    quality_report: pd.DataFrame
    correlation_report: pd.DataFrame


def clean_doi_value(value: object) -> str:
    """Return a normalised DOI string suitable for metadata joins."""

    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    lowered = text.lower()
    lowered = _HTTP_DOI_PREFIX_PATTERN.sub("", lowered)
    lowered = _DOI_PREFIX_PATTERN.sub("", lowered)
    return lowered.strip()


def _normalise_chembl_record(item: Mapping[str, object]) -> dict[str, object]:
    """Extract canonical columns from a raw ChEMBL document payload."""

    journal = item.get("journal_full_title") or item.get("journal")
    return {
        "document_chembl_id": item.get("document_chembl_id"),
        "title": item.get("title"),
        "doi": item.get("doi") or item.get("doi_chembl"),
        "pubmed_id": item.get("pubmed_id"),
        "doc_type": item.get("doc_type"),
        "journal": journal,
        "year": item.get("year"),
    }


def _ensure_absolute_url(url: str, base: str) -> str:
    if url.startswith("http"):
        return url
    return urljoin(base.rstrip("/") + "/", url)


def _update_limit_parameter(url: str, limit: int) -> str:
    parsed = urlsplit(url)
    query = dict(parse_qsl(parsed.query, keep_blank_values=True))
    query["limit"] = str(max(1, limit))
    rebuilt = parsed._replace(query=urlencode(query))
    return urlunsplit(rebuilt)


def _resolve_next_url(meta: object, base: str) -> str:
    if not isinstance(meta, Mapping):
        return ""
    candidate = meta.get("next") or meta.get("next_url")
    if not candidate:
        return ""
    return _ensure_absolute_url(str(candidate), base)


def fetch_chembl_documents(
    limit: int | None,
    *,
    config: Config,
    client: ChemblClient | None = None,
    logger: logging.Logger | None = None,
) -> pd.DataFrame:
    """Fetch raw document rows from ChEMBL respecting ``limit``."""

    if limit is not None and limit < 0:
        msg = "limit must be non-negative"
        raise ValueError(msg)

    chembl_cfg = config.sources.chembl.pipelines.document.chembl
    timeout = chembl_cfg.timeout
    page_size = max(1, min(chembl_cfg.chunk_size, _CHEMBL_PAGE_SIZE_MAX))
    api_cfg = config.api
    base_endpoint = f"{api_cfg.chembl_base.rstrip('/')}/document.json?format=json"

    created_client = False
    if client is None:
        client = ChemblClient(config.api, config.system.retry, config.chembl)
        created_client = True

    try:
        records: list[dict[str, object]] = []
        url = _ensure_absolute_url(f"{base_endpoint}&limit={page_size}", api_cfg.chembl_base)
        while url:
            if limit is not None and len(records) >= limit:
                break
            if limit is not None:
                remaining = limit - len(records)
                if remaining <= 0:
                    break
                if remaining < page_size:
                    url = _update_limit_parameter(url, remaining)
            payload = client.request_json(url, cfg=api_cfg, timeout=timeout)
            items = payload.get("documents") or payload.get("document") or []
            for item in items:
                if not isinstance(item, Mapping):
                    continue
                records.append(_normalise_chembl_record(item))
                if limit is not None and len(records) >= limit:
                    break
            if limit is not None and len(records) >= limit:
                break
            next_url = _resolve_next_url(payload.get("page_meta"), api_cfg.chembl_base)
            if not next_url:
                break
            url = next_url
        frame = pd.DataFrame.from_records(records, columns=_CHEMBL_DOCUMENT_COLUMNS)
        if limit is not None and len(frame) > limit:
            frame = frame.head(limit)
        if logger is not None:
            logger.info(
                "chembl_documents_fetched",
                extra={"rows": len(frame), "limit": limit, "page_size": page_size},
            )
        return frame
    finally:
        if created_client:
            client.close()


def normalize_document_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Return a cleaned frame ready for DOI joins and schema validation."""

    normalized = df.copy(deep=True)
    for column in _CHEMBL_DOCUMENT_COLUMNS:
        if column not in normalized.columns:
            normalized[column] = pd.NA

    def _strip_string(series: pd.Series) -> pd.Series:
        result = series.astype("string").str.strip()
        return result.replace({"": pd.NA})

    normalized["document_chembl_id"] = _strip_string(normalized["document_chembl_id"])
    normalized["title"] = _strip_string(normalized["title"])
    normalized["doc_type"] = _strip_string(normalized["doc_type"])
    normalized["journal"] = _strip_string(normalized["journal"])

    normalized["year"] = (
        pd.to_numeric(normalized["year"], errors="coerce")
        .astype("Int64")
        .astype("string")
    )
    normalized.loc[normalized["year"] == "<NA>", "year"] = pd.NA

    pubmed_numeric = pd.to_numeric(normalized["pubmed_id"], errors="coerce")
    normalized["pubmed_id"] = pubmed_numeric.astype("Int64").astype("string")
    normalized.loc[pubmed_numeric.isna(), "pubmed_id"] = pd.NA

    doi_clean = normalized["doi"].map(clean_doi_value)
    normalized["doi"] = doi_clean.replace({"": pd.NA}).astype("string")
    normalized["doi_key"] = doi_clean.replace({"": pd.NA}).astype("string")

    return normalized


def _prepare_crossref_record(doi: str) -> dict[str, object]:
    return {
        "doi_key": doi,
        "crossref_title": pd.NA,
        "crossref_doc_type": pd.NA,
        "crossref_subject": pd.NA,
        "crossref_error": pd.NA,
    }


def _first_text(entry: object) -> str:
    if isinstance(entry, str):
        return entry
    if isinstance(entry, Sequence) and not isinstance(entry, (str, bytes)):
        for item in entry:
            if isinstance(item, str) and item:
                return item
    return ""


def fetch_crossref_metadata(
    dois: Iterable[str],
    *,
    config: Config,
    logger: logging.Logger | None = None,
    session: Session | None = None,
    limiter: RateLimiter | None = None,
    fetcher: Callable[[str], tuple[Mapping[str, object] | str | None, str]] | None = None,
) -> pd.DataFrame:
    """Retrieve CrossRef metadata for ``dois`` returning a normalised frame."""

    unique = [doi for doi in {clean_doi_value(doi) for doi in dois} if doi]
    if not unique:
        return pd.DataFrame(columns=_CROSSREF_COLUMNS)

    created_session = False
    if session is None:
        session = crossref_session(config.api, config.system.retry, config.crossref)
        created_session = True
    if limiter is None:
        limiter = get_limiter("crossref", config.crossref.rps, config.crossref.burst)

    retry_cfg = config.system.retry

    def _default_fetcher(doi: str) -> tuple[Mapping[str, object] | str | None, str]:
        limiter.acquire()
        data, error = crossref_client.fetch_crossref(
            session,
            doi,
            cfg=config.crossref,
            limiter=limiter,
            retry_cfg=retry_cfg,
        )
        return data, error

    effective_fetcher = fetcher or _default_fetcher

    records: list[dict[str, object]] = []
    try:
        for doi in unique:
            record = _prepare_crossref_record(doi)
            data, error = effective_fetcher(doi)
            if error:
                record["crossref_error"] = error
                records.append(record)
                continue
            if isinstance(data, Mapping):
                message = data.get("message")
                if isinstance(message, Mapping):
                    title = _first_text(message.get("title"))
                    doc_type = message.get("type")
                    subject = message.get("subject")
                    record["crossref_title"] = title or pd.NA
                    record["crossref_doc_type"] = str(doc_type) if doc_type else pd.NA
                    record["crossref_subject"] = (
                        "; ".join(sorted({sub for sub in subject if isinstance(sub, str)}))
                        if isinstance(subject, Sequence)
                        else pd.NA
                    )
                else:
                    record["crossref_error"] = "Invalid response"
            elif isinstance(data, str):
                record["crossref_error"] = data
            else:
                record["crossref_error"] = "Invalid response"
            records.append(record)
    finally:
        if created_session:
            session.close()

    frame = pd.DataFrame.from_records(records, columns=_CROSSREF_COLUMNS)
    for column in frame.columns:
        frame[column] = frame[column].astype("string")
    if logger is not None:
        logger.info("crossref_metadata_fetched", extra={"rows": len(frame)})
    return frame


def _prepare_openalex_record(doi: str | None) -> MutableMapping[str, object]:
    return {
        "doi_key": doi,
        "openalex_title": pd.NA,
        "openalex_doc_type": pd.NA,
        "openalex_type_crossref": pd.NA,
        "openalex_publication_year": pd.NA,
        "openalex_error": pd.NA,
    }


def fetch_openalex_metadata(
    pmids: Iterable[str],
    *,
    config: Config,
    logger: logging.Logger | None = None,
    session: Session | None = None,
    limiter: RateLimiter | None = None,
    fetcher: Callable[[str], tuple[Mapping[str, object] | str | None, str]] | None = None,
) -> pd.DataFrame:
    """Retrieve OpenAlex metadata for ``pmids`` returning DOI-keyed rows."""

    unique_pmids = [pmid for pmid in {str(pmid).strip() for pmid in pmids} if pmid]
    if not unique_pmids:
        return pd.DataFrame(columns=_OPENALEX_COLUMNS)

    created_session = False
    if session is None:
        session = openalex_session(config.api, config.system.retry, config.openalex)
        created_session = True
    if limiter is None:
        limiter = get_limiter("openalex", config.openalex.rps, config.openalex.burst)

    retry_cfg = config.system.retry

    def _default_fetcher(pmid: str) -> tuple[Mapping[str, object] | str | None, str]:
        limiter.acquire()
        data, error = openalex_client.fetch_openalex(
            session,
            pmid,
            cfg=config.openalex,
            limiter=limiter,
            retry_cfg=retry_cfg,
        )
        return data, error

    effective_fetcher = fetcher or _default_fetcher

    records: list[dict[str, object]] = []
    try:
        for pmid in unique_pmids:
            data, error = effective_fetcher(pmid)
            if error:
                record = _prepare_openalex_record(None)
                record["openalex_error"] = error
                records.append(record)
                continue
            if not isinstance(data, Mapping):
                record = _prepare_openalex_record(None)
                record["openalex_error"] = "Invalid response"
                records.append(record)
                continue
            doi_value = clean_doi_value(data.get("doi"))
            record = _prepare_openalex_record(doi_value or None)
            display_name = data.get("display_name") or data.get("title")
            record["openalex_title"] = display_name if display_name else pd.NA
            record["openalex_doc_type"] = (
                str(data.get("type")) if data.get("type") else pd.NA
            )
            record["openalex_type_crossref"] = (
                str(data.get("type_crossref")) if data.get("type_crossref") else pd.NA
            )
            publication_year = data.get("publication_year")
            if isinstance(publication_year, int):
                record["openalex_publication_year"] = str(publication_year)
            records.append(record)
    finally:
        if created_session:
            session.close()

    frame = pd.DataFrame.from_records(records, columns=_OPENALEX_COLUMNS)
    for column in frame.columns:
        frame[column] = frame[column].astype("string")
    frame["doi_key"] = frame["doi_key"].replace({"": pd.NA})
    if logger is not None:
        logger.info("openalex_metadata_fetched", extra={"rows": len(frame)})
    return frame.dropna(subset=["doi_key"], how="all")


def merge_document_metadata(
    base_df: pd.DataFrame,
    *,
    crossref_df: pd.DataFrame | None = None,
    openalex_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return ``base_df`` augmented with CrossRef and OpenAlex metadata."""

    merged = base_df.copy(deep=True)
    if "doi_key" not in merged.columns:
        merged["doi_key"] = merged["doi"].astype("string")

    if crossref_df is not None and not crossref_df.empty:
        deduplicated_crossref = crossref_df.drop_duplicates(subset=["doi_key"])
        merged = merged.merge(
            deduplicated_crossref,
            how="left",
            on="doi_key",
        )
    else:
        for column in _CROSSREF_COLUMNS[1:]:
            merged[column] = pd.NA

    if openalex_df is not None and not openalex_df.empty:
        deduplicated_openalex = openalex_df.drop_duplicates(subset=["doi_key"])
        merged = merged.merge(
            deduplicated_openalex,
            how="left",
            on="doi_key",
        )
    else:
        for column in _OPENALEX_COLUMNS[1:]:
            merged[column] = pd.NA

    if "crossref_title" in merged.columns:
        merged["title"] = merged["title"].combine_first(merged["crossref_title"])
    if "openalex_title" in merged.columns:
        merged["title"] = merged["title"].combine_first(merged["openalex_title"])

    if "crossref_doc_type" in merged.columns:
        merged["doc_type"] = merged["doc_type"].combine_first(merged["crossref_doc_type"])
    if "openalex_doc_type" in merged.columns:
        merged["doc_type"] = merged["doc_type"].combine_first(merged["openalex_doc_type"])

    return merged


def _build_quality_artifacts(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    profiler = TableQualityProfiler()
    profiler.consume(frame)
    return build_reports_from_profiler(profiler)


def fetch_normalize_document(
    limit: int | None,
    config: Config,
    logger: logging.Logger,
) -> DocumentFetchResult:
    """Fetch, enrich, and validate ChEMBL documents using DOI metadata."""

    base_df = fetch_chembl_documents(limit, config=config, logger=logger)
    normalized = normalize_document_frame(base_df)

    dois = normalized["doi_key"].dropna().unique().tolist()
    crossref_df = fetch_crossref_metadata(dois, config=config, logger=logger)

    pmids = normalized["pubmed_id"].dropna().unique().tolist()
    openalex_df = fetch_openalex_metadata(pmids, config=config, logger=logger)

    merged = merge_document_metadata(
        normalized,
        crossref_df=crossref_df,
        openalex_df=openalex_df,
    )

    validated = validate_document_frame(merged)

    quality_report, correlation_report = _build_quality_artifacts(validated)

    return DocumentFetchResult(
        data=validated,
        quality_report=quality_report,
        correlation_report=correlation_report,
    )


__all__ = [
    "DocumentFetchResult",
    "clean_doi_value",
    "fetch_chembl_documents",
    "fetch_crossref_metadata",
    "fetch_normalize_document",
    "fetch_openalex_metadata",
    "merge_document_metadata",
    "normalize_document_frame",
]
