from __future__ import annotations

from pathlib import Path
import re

import pytest


_MERGE_CONFLICT_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^<<<<<<< ", re.MULTILINE),
    re.compile(r"^>>>>>>> ", re.MULTILINE),
)


def _iter_python_sources(root: Path) -> list[Path]:
    """Return tracked Python sources excluding transient caches."""

    sources: list[Path] = []
    for path in root.rglob("*.py"):
        relative = path.relative_to(root)
        parts = relative.parts
        if any(part.startswith(".") for part in parts):
            continue
        if any(part == "__pycache__" for part in parts):
            continue
        sources.append(path)
    return sources


@pytest.mark.unit
def test_python_sources__have_no_merge_conflict_markers() -> None:
    root = Path(__file__).resolve().parents[2]
    offenders: list[str] = []

    for source in _iter_python_sources(root):
        text = source.read_text(encoding="utf-8")
        if any(pattern.search(text) for pattern in _MERGE_CONFLICT_PATTERNS):
            offenders.append(source.relative_to(root).as_posix())

    assert not offenders, (
        "Merge conflict markers detected in tracked sources: "
        + ", ".join(sorted(offenders))
    )
