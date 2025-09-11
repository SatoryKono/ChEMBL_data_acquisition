"""Pydantic models for CSV export CLI arguments.

This module defines :class:`CSVExportArgs`, a Pydantic model that validates
and normalises arguments used by :mod:`csv_utils_main`. Using a data model
ensures consistent handling of command-line parameters across scripts.
"""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field


class CSVExportArgs(BaseModel):  # type: ignore[misc]
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
    log_level:
        Logging verbosity passed through to the application logger.
    """

    model_config = ConfigDict(extra="forbid")

    input_csv: Path = Field(default=Path("input.csv"))
    output_csv: Path | None = None
    sep: str = ","
    encoding: str = "utf8"
    col_order: list[str] | None = None
    key_cols: list[str] | None = None
    log_level: str = "INFO"


__all__ = ["CSVExportArgs"]
