"""Lightweight Git helpers.

This module provides shared utilities related to Git.
Currently it exposes :func:`_git_sha` which returns the current
commit hash.  The helper is centralised here so other modules can
record provenance without duplicating logic.
"""

from __future__ import annotations

import functools
import os
import shutil
import subprocess
from pathlib import Path

from .log import logger


@functools.lru_cache(maxsize=1)
def _git_sha() -> str:
    """Return the Git commit hash for the repository.

    The command is executed only once and the result cached to avoid
    repeated subprocess calls.  If Git is unavailable or the repository
    lacks its ``.git`` directory, ``"UNKNOWN"`` is returned and a warning
    is logged.

    Returns
    -------
    str
        The commit hash or ``"UNKNOWN"`` when it cannot be determined.
    """

    repo_root = Path(__file__).resolve().parent.parent
    env_sha = os.getenv("GIT_SHA")
    if env_sha:
        logger.info("git_sha_env", sha=env_sha)
        return env_sha

    git_executable = shutil.which("git")
    if git_executable is None:
        logger.warning("git_executable_missing")
        return "UNKNOWN"

    if not (repo_root / ".git").exists():
        logger.warning("git_directory_missing", path=str(repo_root))
        return "UNKNOWN"

    try:
        result = subprocess.check_output(
            [git_executable, "rev-parse", "HEAD"], cwd=repo_root
        )
        return result.decode().strip()
    except subprocess.CalledProcessError as exc:  # pragma: no cover - rare
        logger.warning("git_sha_unavailable", extra={"error": str(exc)})
        return "UNKNOWN"
