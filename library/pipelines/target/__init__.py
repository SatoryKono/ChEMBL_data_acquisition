"""Target acquisition pipeline helpers and orchestrators."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

from ...config import Config
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

__all__ = [
    "TargetPipelineOptions",
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


@dataclass(slots=True)
class TargetPipelineOptions:
    """Configuration for executing the combined target pipeline."""

    input_csv: Path
    output_csv: Path
    command: str = "all"
    limit: int | None = None
    offset: int = 0
    skip_existing: bool = False
    force: bool = False
    chunk_size: int | None = None
    timeout: float | None = None


def run_pipeline(config: Config, options: TargetPipelineOptions) -> int:
    """Execute the requested target pipeline mode using CLI helpers."""

    from scripts import get_target_data as _target_cli

    command = (options.command or "all").strip().lower()
    handler_map = {
        "chembl": _target_cli.run_chembl,
        "uniprot": _target_cli.run_uniprot,
        "iuphar": _target_cli.run_iuphar,
        "all": _target_cli.run_all,
    }
    handler = handler_map.get(command)
    if handler is None:
        raise ValueError(f"unsupported target pipeline command: {command}")

    args = argparse.Namespace(
        input_csv=Path(options.input_csv),
        output_csv=Path(options.output_csv),
        final_out=Path(options.output_csv),
        limit=options.limit,
        offset=options.offset,
        skip_existing=options.skip_existing,
        force=options.force,
        command=command,
        chunk_size=options.chunk_size,
        timeout=options.timeout,
        func=handler,
    )

    if options.skip_existing and options.output_csv.exists() and not options.force:
        return 0

    return int(handler(config, args))
