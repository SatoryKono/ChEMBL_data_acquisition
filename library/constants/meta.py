"""Schema for CSV metadata sidecar files."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, Field


class CsvMetaSchema(BaseModel):
    """Pydantic model describing ``.meta.yaml`` sidecars.

    Attributes
    ----------
    generated_at:
        Timestamp when the export was created.
    git_sha:
        Git commit hash of the repository at creation time.
    columns:
        Ordered list of column names present in the CSV file.
    dtypes:
        Mapping of column names to their dtypes as reported by pandas.
    command:
        Optional command line invocation string.
    config:
        Optional serialised configuration used to generate the file.
    """

    generated_at: datetime
    git_sha: str
    columns: list[str] = Field(default_factory=list)
    dtypes: dict[str, str] = Field(default_factory=dict)
    command: str | None = None
    config: dict[str, Any] = Field(default_factory=dict)


__all__ = ["CsvMetaSchema"]
