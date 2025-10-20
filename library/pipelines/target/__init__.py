"""Target acquisition pipeline helpers and programmatic runner."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
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
    date: str | None = None
    output_stem: str | None = None


def run_pipeline(config: Config, options: TargetPipelineOptions) -> PipelineRunResult:
    """Execute the target pipeline for the selected command."""

    from library.cli.commands import get_target_data as target_cli

    return target_cli.run_target_service(config, options)


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
