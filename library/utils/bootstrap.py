"""Runtime helpers for script execution.

This module ensures that executing a script via ``python scripts/foo.py``
behaves the same as ``python -m scripts.foo`` by pre-pending the project
root to :mod:`sys.path`.
"""

from __future__ import annotations

import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]


def ensure_project_root() -> None:
  """Insert the project root into :data:`sys.path` if missing."""
  project_root = str(_PROJECT_ROOT)
  if project_root not in sys.path:
    sys.path.insert(0, project_root)
