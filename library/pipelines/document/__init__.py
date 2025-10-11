"""Document aggregation and enrichment pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from library.config import Config
from library.pipelines.common import PipelineRunResult

from . import pipeline, postprocessing, type_classifier, type_terms
from .service import run_document_service
from .chembl_document import get_documents


@dataclass(slots=True)
class DocumentPipelineOptions:
    """Configuration for executing the document pipeline programmatically.

    The :mod:`library.cli.commands.get_data` orchestration forwards ``skip_existing``
    when the ``--skip-existing`` flag is set so repeated runs can reuse previously
    generated CSV files instead of rewriting them.
    """

    input_csv: Path
    output_csv: Path
    mode: Literal["chembl", "pubmed", "all"] = "all"
    limit: int | None = None
    offset: int = 0
    force: bool = False
    skip_existing: bool = False
    date_prefix: str | None = None
    fallback_doi_enabled: bool = False
    fallback_doi_path: Path | None = None
    fallback_doi_overwrite: bool = False
    fallback_doi_delimiter: str | None = None
    fallback_doi_encoding: str | None = None
    fallback_doi_col_pmid: str = "PMID"
    fallback_doi_col_doi: str = "DOI"
    rerun_postprocess: bool = False
    date_prefix: str | None = None
    output_stem: str | None = None

def run_pipeline(config: Config, options: DocumentPipelineOptions) -> PipelineRunResult:
    """Run the document pipeline using programmatic options."""

    return run_document_service(config, options)


__all__ = [
    "DocumentPipelineOptions",
    "PipelineRunResult",
    "get_documents",
    "pipeline",
    "postprocessing",
    "run_pipeline",
    "type_classifier",
    "type_terms",
]
