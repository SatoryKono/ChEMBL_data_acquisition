"""Utilities for aggregating validation errors into a sidecar CSV file."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


class SidecarErrors:
    """Collect tabular error records and persist them to CSV.

    The class stores each error as a mapping of column names to values.
    Errors are written to disk only when :meth:`save` is called and at least
    one error was recorded.
    """

    def __init__(self) -> None:
        """Initialize an empty error collection."""
        self._errors: list[dict[str, Any]] = []

    def add_error(self, row: dict[str, Any]) -> None:
        """Add a validation error description.

        Parameters
        ----------
        row: dict[str, Any]
            Mapping describing a single validation failure.
        """
        self._errors.append(row)

    def save(self, path: Path) -> None:
        """Write collected errors to ``path`` as CSV.

        Parameters
        ----------
        path: pathlib.Path
            Destination file. Parent directories are created as needed.

        Notes
        -----
        The file is created only if at least one error was recorded.
        """
        if not self._errors:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = sorted({k for row in self._errors for k in row.keys()})
        with path.open("w", newline="", encoding="utf8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self._errors)
