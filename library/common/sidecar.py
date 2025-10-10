"""Utilities for aggregating validation errors into a sidecar CSV file."""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path
from typing import Any

from ..config import Config
from ..io.metadata import write_meta_yaml


class SidecarErrors:
    """Collect tabular error records and persist them to CSV.

    The class stores each error as a mapping of column names to values.
    Errors are written to disk only when :meth:`save` is called and at least
    one error was recorded.  When an optional ``chunk_size`` is provided the
    class will periodically spill buffered rows to a temporary CSV file in
    order to bound memory usage.
    """

    def __init__(self, *, chunk_size: int | None = None) -> None:
        """Initialize an empty error collection.

        Parameters
        ----------
        chunk_size:
            Optional maximum number of error rows kept in memory.  When the
            buffer reaches this size it is flushed to an on-disk CSV.
        """
        if chunk_size is not None and chunk_size <= 0:
            msg = "chunk_size must be a positive integer"
            raise ValueError(msg)
        self._chunk_size = chunk_size
        self._errors: list[dict[str, Any]] = []
        self._fieldnames: set[str] = set()
        self._overflow_path: Path | None = None
        self._headers_written = False
        self._needs_rewrite = False

    def add_error(self, row: dict[str, Any]) -> None:
        """Add a validation error description.

        Parameters
        ----------
        row: dict[str, Any]
            Mapping describing a single validation failure.
        """
        self._errors.append(row)
        new_fields = set(row.keys()) - self._fieldnames
        if new_fields:
            self._fieldnames.update(new_fields)
            if self._headers_written:
                self._needs_rewrite = True
        if self._chunk_size is not None and len(self._errors) >= self._chunk_size:
            self._flush_buffer()

    def save(self, path: Path, *, cfg: Config | None = None) -> None:
        """Write collected errors to ``path`` as CSV and emit metadata.

        Parameters
        ----------
        path:
            Destination file. Parent directories are created as needed.
        cfg:
            Optional configuration forwarded to :func:`write_meta_yaml`.

        Notes
        -----
        The file is created only if at least one error was recorded.
        """
        if not self._errors and self._overflow_path is None:
            return

        path.parent.mkdir(parents=True, exist_ok=True)

        if self._overflow_path is not None and self._errors:
            self._flush_buffer()

        if not self._fieldnames:
            self._fieldnames.update({k for row in self._errors for k in row.keys()})

        fieldnames = sorted(self._fieldnames)

        if self._overflow_path is None:
            with path.open("w", newline="", encoding="utf8") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, restval="")
                writer.writeheader()
                writer.writerows(self._errors)
        else:
            if self._needs_rewrite:
                self._rewrite_overflow()
            with path.open("w", newline="", encoding="utf8") as dst:
                writer = csv.DictWriter(dst, fieldnames=fieldnames, restval="")
                writer.writeheader()
                with self._overflow_path.open("r", newline="", encoding="utf8") as src:
                    reader = csv.DictReader(src)
                    for row in reader:
                        writer.writerow(row)
            self._overflow_path.unlink()
            self._overflow_path = None
            self._headers_written = False
            self._needs_rewrite = False

        write_meta_yaml(path, cfg=cfg, columns=fieldnames)
        self._errors.clear()
        self._fieldnames.clear()

    def _flush_buffer(self) -> None:
        if not self._errors:
            return
        if self._overflow_path is None:
            fd, name = tempfile.mkstemp(prefix="sidecar_", suffix=".csv")
            os.close(fd)
            self._overflow_path = Path(name)
        if self._needs_rewrite:
            self._rewrite_overflow()
        fieldnames = sorted(self._fieldnames)
        mode = "a" if self._headers_written else "w"
        with self._overflow_path.open(mode, newline="", encoding="utf8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, restval="")
            if not self._headers_written:
                writer.writeheader()
                self._headers_written = True
            writer.writerows(self._errors)
        self._errors.clear()

    def _rewrite_overflow(self) -> None:
        assert self._overflow_path is not None
        fd, name = tempfile.mkstemp(prefix="sidecar_rewrite_", suffix=".csv")
        os.close(fd)
        tmp_path = Path(name)
        fieldnames = sorted(self._fieldnames)
        with (
            self._overflow_path.open("r", newline="", encoding="utf8") as src,
            tmp_path.open("w", newline="", encoding="utf8") as dst,
        ):
            reader = csv.DictReader(src)
            writer = csv.DictWriter(dst, fieldnames=fieldnames, restval="")
            writer.writeheader()
            for row in reader:
                writer.writerow(row)
        self._overflow_path.unlink()
        tmp_path.rename(self._overflow_path)
        self._needs_rewrite = False
