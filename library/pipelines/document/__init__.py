"""Document aggregation and enrichment pipelines."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Literal

from library.config import Config
from library.pipelines.common import PipelineRunResult

from . import pipeline, postprocessing, type_classifier, type_terms
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
    fallback_doi_enabled: bool = False
    fallback_doi_path: Path | None = None
    fallback_doi_overwrite: bool = False
    fallback_doi_delimiter: str | None = None
    fallback_doi_encoding: str | None = None
    fallback_doi_col_pmid: str = "PMID"
    fallback_doi_col_doi: str = "DOI"


def _update_document_section(
    cfg: Config,
    mode: Literal["chembl", "pubmed", "all"],
    *,
    limit: int | None,
    offset: int,
) -> None:
    section = getattr(cfg.sources.chembl.pipelines.document, mode)
    updates: dict[str, object] = {"offset": offset}
    if limit is not None:
        updates["limit"] = limit
    setattr(
        cfg.sources.chembl.pipelines.document,
        mode,
        section.model_copy(update=updates),
    )


def run_pipeline(config: Config, options: DocumentPipelineOptions) -> PipelineRunResult:
    """Run the document pipeline using programmatic options."""

    output_path = Path(options.output_csv)
    if options.skip_existing and output_path.exists() and not options.force:
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=False,
            reason="skip_existing",
            written=False,
        )

    from scripts import get_document_data as document_cli  # Lazy import to avoid cycles

    cfg = config.model_copy(deep=True)
    _update_document_section(
        cfg,
        options.mode,
        limit=options.limit,
        offset=options.offset,
    )

    if options.mode == "chembl":
        runner = document_cli.run_chembl
    elif options.mode == "pubmed":
        runner = document_cli.run_pubmed
    elif options.mode == "all":
        runner = document_cli.run_all
    else:  # pragma: no cover - defensive guard
        msg = f"unsupported document pipeline mode: {options.mode}"
        raise ValueError(msg)

    args = SimpleNamespace(
        input_csv=Path(options.input_csv),
        final_out=output_path,
        output_csv=output_path,
        skip_existing=options.skip_existing,
        force=options.force,
        limit=options.limit,
        offset=options.offset,
        fallback_doi_enabled=options.fallback_doi_enabled,
        fallback_doi_path=options.fallback_doi_path,
        fallback_doi_overwrite=options.fallback_doi_overwrite,
        fallback_doi_delimiter=options.fallback_doi_delimiter,
        fallback_doi_encoding=options.fallback_doi_encoding,
        fallback_doi_col_pmid=options.fallback_doi_col_pmid,
        fallback_doi_col_doi=options.fallback_doi_col_doi,
        mode=options.mode,
        command=options.mode,
    )

    exit_code = int(runner(cfg, args))
    reason = None if exit_code == 0 else "pipeline_failed"
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=exit_code == 0,
    )


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
