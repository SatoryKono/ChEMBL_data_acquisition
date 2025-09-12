"""Helpers for writing metadata sidecar files.

This module provides :func:`write_meta_yaml` which records execution context
and output statistics in a YAML file alongside a generated CSV.  The metadata
captures information such as the Git commit, Python version, and dataset
summary statistics.
"""

from __future__ import annotations

import functools
import hashlib
import os
import platform
import shutil
import subprocess
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypedDict

import yaml

from .config import _mask_secrets
from .log import logger

# ``datetime.UTC`` is only available in Python 3.11 and later.
# ``timezone.utc`` provides the same value and works on older versions.
UTC = timezone.utc  # noqa: UP017


class Stats(TypedDict):
    """Execution statistics written to the metadata file."""

    rows_total: int
    rows_kept: int
    rows_dropped: int
    output_sha256: str


def file_sha256(path: Path | str) -> str:
    """Return the SHA256 checksum of ``path``.

    Parameters
    ----------
    path:
        File path for which the digest should be computed.

    Returns
    -------
    str
        Hex encoded SHA256 digest of the file contents.
    """

    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@functools.lru_cache(maxsize=1)
def _git_sha() -> str:
    """Return the current Git commit hash.


    The result is cached to avoid repeated invocations of Git. If Git is
    unavailable, ``"UNKNOWN"`` is returned and a warning is logged.

    """

    repo_root = Path(__file__).resolve().parent.parent
    env_sha = os.getenv("GIT_SHA")
    if env_sha:
        logger.info("Using git SHA from GIT_SHA: %s", env_sha)
        return env_sha

    try:
        result = subprocess.check_output(
            [git_executable, "rev-parse", "HEAD"], cwd=repo_root
        )
        return result.decode().strip()
    except subprocess.CalledProcessError as exc:
        logger.warning("Unable to determine git SHA: %s", exc)
        return "UNKNOWN"


def write_meta_yaml(
    csv_path: Path | str,
    command: str,
    config_subset: Mapping[str, Any],
    inputs: Mapping[str, Any],
    stats: Stats,
    schema: str,
) -> Path:
    """Write metadata for ``csv_path`` to ``<csv_path>.meta.yaml``.

    Parameters
    ----------
    csv_path:
        Path to the generated CSV file.
    command:
        Command used to invoke the pipeline or script.
    config_subset:
        Relevant configuration values. Secret keys are masked before
        serialisation.
    inputs:
        Description of the inputs that produced the output.
    stats:
        Summary statistics about the output table.
    schema:
        Name of the schema applied to the output data.

    Returns
    -------
    Path
        Path to the written metadata YAML file.
    """

    path = Path(csv_path)
    meta_path = path.with_name(path.name + ".meta.yaml")

    metadata: dict[str, Any] = {
        "generated_at": datetime.now(UTC).isoformat(),
        "git_sha": _git_sha(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "command": command,
        "config": _mask_secrets(dict(config_subset)),
        "inputs": dict(inputs),
        "stats": stats,
        "schema": schema,
    }

    with meta_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(metadata, fh, sort_keys=False)
        logger.info("Metadata written to %s", meta_path)

    return meta_path
