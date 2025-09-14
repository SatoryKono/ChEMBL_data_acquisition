"""Pydantic models for CSV export CLI arguments.

This module defines :class:`CSVExportArgs`, a Pydantic model that validates
and normalises arguments used by :mod:`scripts.csv_utils_main`. Using a data model
ensures consistent handling of command-line parameters across scripts.
"""

# ruff: noqa: UP007,UP045

from __future__ import annotations

from pathlib import Path
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class CSVExportArgs(BaseModel):
    """Validated parameters for CSV export operations.

    Attributes
    ----------
    input_csv:
        Path to the source CSV file.
    output_csv:
        Optional destination path. When ``None`` an automatic name based on
        ``input_csv`` is used.
    sep:
        Field delimiter for the CSV file.
    encoding:
        Character encoding of the CSV file.
    col_order:
        Optional list specifying the expected column order. Columns not listed
        here are dropped during export.
    key_cols:
        Optional list of columns used to determine row ordering.
    chunk_size:
        Number of rows processed per chunk when streaming CSV input.
    log_level:
        Logging verbosity passed through to the application logger.
    """

    model_config = ConfigDict(extra="forbid")

    input_csv: Path = Field(default=Path("input.csv"))
    output_csv: Optional[Path] = None  # noqa: UP007
    sep: str = ","
    encoding: str = "utf8"
    col_order: Optional[list[str]] = None  # noqa: UP007
    key_cols: Optional[list[str]] = None  # noqa: UP007
    chunk_size: int = 1000
    log_level: str = "INFO"


__all__ = ["CSVExportArgs"]
