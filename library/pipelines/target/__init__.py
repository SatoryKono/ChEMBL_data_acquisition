"""Target acquisition pipeline helpers and programmatic runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Literal

from library.config import Config
from library.pipelines.common import PipelineRunResult

from . import (
    cellularity,
    helpers,
    multifunctional,
    organism_classification,
    pipeline,
    postprocessing,
    protein_classification,
)
from .chembl_target import normalize_reaction_ec_numbers


@dataclass(slots=True)
class TargetPipelineOptions:
    """Configuration options for executing the target pipeline.

    ``skip_existing`` is set by the unified ``get_data`` CLI when ``--skip-existing``
    is requested so cached target exports are preserved during incremental runs.
    """

    input_csv: Path
    output_csv: Path
    command: Literal["chembl", "uniprot", "iuphar", "all"] = "all"
    limit: int | None = None
    offset: int = 0
    raw_output: Path | None = None
    raw_format: str = "csv"
    no_reindex_raw: bool = False
    id_columns: tuple[str, ...] | None = None
    skip_existing: bool = False
    force: bool = False


def _update_target_config(
    cfg: Config,
    *,
    command: str,
    limit: int | None,
    offset: int,
) -> None:
    target_cfg = cfg.sources.chembl.pipelines.target

    def _apply(section_name: str) -> None:
        section = getattr(target_cfg, section_name)
        updates: dict[str, object] = {"offset": offset}
        if limit is not None:
            updates["limit"] = limit
        setattr(target_cfg, section_name, section.model_copy(update=updates))

    if command == "all":
        _apply("all")
        _apply("chembl")
        _apply("uniprot")
        _apply("iuphar")
    else:
        _apply(command)


def run_pipeline(config: Config, options: TargetPipelineOptions) -> PipelineRunResult:
    """Execute the target pipeline for the selected command."""

    output_path = Path(options.output_csv)
    if options.skip_existing and output_path.exists() and not options.force:
        return PipelineRunResult(
            exit_code=0,
            output_path=output_path,
            executed=False,
            reason="skip_existing",
            written=False,
        )

    from scripts import get_target_data as target_cli  # Lazy import to avoid cycles

    cfg = config.model_copy(deep=True)
    _update_target_config(
        cfg,
        command=options.command,
        limit=options.limit,
        offset=options.offset,
    )

    args = SimpleNamespace(
        input_csv=Path(options.input_csv),
        final_out=output_path,
        output_csv=output_path,
        raw_out=options.raw_output,
        raw_format=options.raw_format,
        no_reindex_raw=options.no_reindex_raw,
        id_cols=options.id_columns,
        limit=options.limit,
        offset=options.offset,
        command=options.command,
        skip_existing=False,
        force=options.force,
    )

    if options.command == "chembl":
        runner = target_cli.run_chembl
    elif options.command == "uniprot":
        runner = target_cli.run_uniprot
    elif options.command == "iuphar":
        runner = target_cli.run_iuphar
    elif options.command == "all":
        runner = target_cli.run_all
    else:  # pragma: no cover - defensive guard
        msg = f"unsupported target pipeline command: {options.command}"
        raise ValueError(msg)

    exit_code = int(runner(cfg, args))
    reason = None if exit_code == 0 else "pipeline_failed"
    written = None if exit_code != 0 else True
    return PipelineRunResult(
        exit_code=exit_code,
        output_path=output_path,
        executed=True,
        reason=reason,
        written=written,
    )


__all__ = [
    "TargetPipelineOptions",
    "PipelineRunResult",
    "normalize_reaction_ec_numbers",
    "cellularity",
    "helpers",
    "multifunctional",
    "organism_classification",
    "pipeline",
    "postprocessing",
    "protein_classification",
    "run_pipeline",
]
