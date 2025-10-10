"""Data structures shared by pipeline execution helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class PipelineRunResult:
    """Summarise the outcome of a programmatic pipeline invocation."""

    exit_code: int
    output_path: Path
    executed: bool = True
    reason: str | None = None
    written: bool | None = None
