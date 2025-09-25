from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def smoke_output_dir() -> Path:
    """Return the directory for smoke-test outputs."""

    path = Path("data/output-smoke")
    path.mkdir(parents=True, exist_ok=True)
    return path
