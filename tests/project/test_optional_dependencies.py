"""Tests for optional dependency configuration.

Ensures each optional dependency list in the project metadata remains
free of duplicate requirement entries.
"""

from __future__ import annotations

import tomllib
from collections import defaultdict
from pathlib import Path


def _project_root() -> Path:
    path = Path(__file__).resolve().parent
    while path.name != "tests":
        path = path.parent
    return path.parent


def test_optional_dependencies_are_unique() -> None:
    """Verify that optional dependency groups contain unique entries."""

    project_root = _project_root()
    pyproject_path = project_root / "pyproject.toml"
    data = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))

    duplicates: dict[str, set[str]] = defaultdict(set)

    for extra_name, requirements in (
        data.get("project", {}).get("optional-dependencies", {}).items()
    ):
        seen: set[str] = set()
        for requirement in requirements:
            normalized = requirement.strip()
            if normalized in seen:
                duplicates[extra_name].add(normalized)
            else:
                seen.add(normalized)

    assert (
        not duplicates
    ), "Duplicate optional dependency entries detected: " + ", ".join(
        f"{group}: {sorted(values)}" for group, values in sorted(duplicates.items())
    )
