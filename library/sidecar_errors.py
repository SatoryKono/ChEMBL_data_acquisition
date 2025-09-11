"""Utilities for writing validation errors to sidecar files."""

from __future__ import annotations

from pathlib import Path
from typing import List


class SidecarErrors:
    """Collect and persist error messages to a sidecar file.

    Parameters
    ----------
    path:
        Location of the sidecar file where error messages will be written.
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self._errors: List[str] = []

    def add(self, message: str) -> None:
        """Record an error ``message`` for later persistence."""
        self._errors.append(message)

    def write(self) -> None:
        """Write collected error messages to ``self.path``.

        The file is created only if at least one message was recorded.
        """
        if not self._errors:
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text("\n".join(self._errors) + "\n", encoding="utf8")
