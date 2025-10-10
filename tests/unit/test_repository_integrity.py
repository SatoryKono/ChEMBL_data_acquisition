from __future__ import annotations
import py_compile
import re
import tempfile
from pathlib import Path

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


@pytest.mark.unit
def test_python_sources__compile_without_syntax_errors() -> None:
    """Ensure all tracked Python modules compile successfully."""

    root = Path(__file__).resolve().parents[2]
    failures: list[str] = []

    with tempfile.TemporaryDirectory(prefix="chembl_compile_") as tmpdir:
        cache_dir = Path(tmpdir)
        for index, source in enumerate(_iter_python_sources(root)):
            cfile = cache_dir / f"module_{index}.pyc"
            try:
                py_compile.compile(
                    str(source),
                    cfile=str(cfile),
                    doraise=True,
                )
            except py_compile.PyCompileError as exc:  # pragma: no cover - exercised on failure
                failures.append(
                    f"{source.relative_to(root).as_posix()}: {exc.exc_type_name}: {exc.msg}"
                )
            finally:
                if cfile.exists():
                    try:
                        cfile.unlink()
                    except OSError:
                        pass

    assert not failures, (
        "Python sources with syntax errors detected: "
        + ", ".join(sorted(failures))
    )
