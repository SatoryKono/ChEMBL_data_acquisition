"""Helpers for writing metadata sidecar files.

This module provides :func:`write_meta_yaml` which records execution context
and output statistics in a YAML file alongside a generated CSV.  The metadata
captures information such as the Git commit, Python version, and dataset
summary statistics.
"""

from __future__ import annotations

import hashlib
import platform
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypedDict

import yaml

from .config import _mask_secrets
from .git_utils import _git_sha
from .log import logger

# ``datetime.UTC`` is only available in Python 3.11 and later.
# ``timezone.utc`` provides the same value and works on older versions.
UTC = timezone.utc  # noqa: UP017


class _StatsRequired(TypedDict):
    rows_total: int
    rows_kept: int
    rows_dropped: int
    output_sha256: str


class Stats(_StatsRequired, total=False):
    """Execution statistics written to the metadata file."""

    parent_lookup_source: str
    parent_lookup_missing: int
    missing_molecule_ids: list[str]
    missing_molecule_ids_count: int


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

    existing: dict[str, Any] = {}
    if meta_path.exists():
        try:
            with meta_path.open("r", encoding="utf-8") as fh:
                loaded = yaml.safe_load(fh)
        except (OSError, yaml.YAMLError):
            loaded = None
        if isinstance(loaded, dict):
            existing = dict(loaded)

    metadata: dict[str, Any] = dict(existing)
    metadata.update(
        {
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
    )

    with meta_path.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(metadata, fh, sort_keys=False)
        logger.info("metadata_written", path=str(meta_path))

    return meta_path
