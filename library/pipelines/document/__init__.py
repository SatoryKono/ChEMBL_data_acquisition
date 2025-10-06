"""Document aggregation and enrichment pipelines."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from ...config import Config
from . import pipeline, postprocessing, type_classifier, type_terms
from .chembl_document import get_documents

__all__ = [
    "DocumentPipelineOptions",
    "get_documents",
    "pipeline",
    "postprocessing",
    "run_pipeline",
    "type_classifier",
    "type_terms",
]


@dataclass(slots=True)
class DocumentPipelineOptions:
    """Settings controlling the document pipeline execution."""

    input_csv: Path
    output_csv: Path
    mode: str = "all"
    limit: int | None = None
    offset: int = 0
    skip_existing: bool = False
    force: bool = False
    timeout: float | None = None
    chembl_timeout: float | None = None
    chembl_chunk_size: int | None = None
    pubmed_workers: int | None = None
    pubmed_batch_size: int | None = None
    pubmed_sleep: float | None = None
    fallback_doi_enabled: bool = False
    fallback_doi_path: Path | None = None
    fallback_doi_overwrite: bool = False
    fallback_doi_delimiter: str | None = None
    fallback_doi_encoding: str | None = None
    fallback_doi_col_pmid: str | None = None
    fallback_doi_col_doi: str | None = None
    column: str | None = None


def run_pipeline(config: Config, options: DocumentPipelineOptions) -> int:
    """Execute the requested document pipeline mode using the CLI helpers."""

    from scripts import get_document_data as _document_cli

    mode = (options.mode or "all").strip().lower()
    handler = _document_cli.MODE_HANDLERS.get(mode)
    if handler is None:
        raise ValueError(f"unsupported document pipeline mode: {mode}")

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        output_csv=Path(options.output_csv),
        final_out=Path(options.output_csv),
        limit=options.limit,
        offset=options.offset,
        skip_existing=options.skip_existing,
        force=options.force,
        mode=mode,
        command=mode,
        timeout=options.timeout,
        chembl_timeout=options.chembl_timeout,
        chembl_chunk_size=options.chembl_chunk_size,
        pubmed_workers=options.pubmed_workers,
        pubmed_batch_size=options.pubmed_batch_size,
        pubmed_sleep=options.pubmed_sleep,
        fallback_doi_enabled=options.fallback_doi_enabled,
        fallback_doi_path=options.fallback_doi_path,
        fallback_doi_overwrite=options.fallback_doi_overwrite,
        fallback_doi_delimiter=options.fallback_doi_delimiter,
        fallback_doi_encoding=options.fallback_doi_encoding,
        fallback_doi_col_pmid=options.fallback_doi_col_pmid,
        fallback_doi_col_doi=options.fallback_doi_col_doi,
        column=options.column,
    )

    timeout_value: float | None = None
    if mode == "chembl":
        timeout_value = options.timeout
    elif mode == "all":
        timeout_value = options.chembl_timeout or options.timeout
    if timeout_value is not None:
        config.api.timeout_read = timeout_value

    if options.skip_existing and options.output_csv.exists() and not options.force:
        return 0

    return int(handler(config, args))
