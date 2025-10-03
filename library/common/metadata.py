"""Helpers for writing metadata sidecar files.

This module provides :func:`write_meta_yaml` which records execution context
and output statistics in a YAML file alongside a generated CSV.  The metadata
captures information such as the Git commit, Python version, and dataset
summary statistics.
"""

from __future__ import annotations

import hashlib
import platform
from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypedDict

import yaml

from ..config import _mask_secrets
from .git import _git_sha
from .log import logger
from ..utils.atomic import open_atomic


def _load_metadata(meta_path: Path) -> dict[str, Any]:
    """Return the metadata dictionary stored in ``meta_path`` if available."""

    if not meta_path.exists():
        return {}

    try:
        with meta_path.open("r", encoding="utf-8") as fh:
            loaded = yaml.safe_load(fh)
    except (OSError, yaml.YAMLError):
        loaded = None

    if isinstance(loaded, dict):
        return dict(loaded)
    return {}

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
    parent_lookup_hierarchy_attached: int
    parent_lookup_fallback_attached: int
    parent_lookup_no_parent: int
    missing_molecule_ids: list[str]
    missing_molecule_ids_count: int
    chunk_fetch_failure_chunks: int
    chunk_fetch_failure_ids: list[str]


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
    *,
    invocation: Sequence[str] | None = None,
    extra_metadata: Mapping[str, Any] | None = None,
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
    invocation:
        Optional tuple describing the exact CLI invocation. Persisted when
        provided to aid reproducibility.
    extra_metadata:
        Optional mapping merged into the generated metadata prior to
        serialisation.

    Returns
    -------
    Path
        Path to the written metadata YAML file.
    """

    path = Path(csv_path)
    meta_path = path.with_name(path.name + ".meta.yaml")

    existing = _load_metadata(meta_path)

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
    if invocation:
        metadata["invocation"] = list(invocation)
    if extra_metadata:
        metadata.update(dict(extra_metadata))

    with open_atomic(meta_path, encoding="utf-8") as fh:
        yaml.safe_dump(metadata, fh, sort_keys=False)
        logger.info("metadata_written", path=str(meta_path))

    return meta_path


def record_quality_failure(
    meta_path: Path | str,
    *,
    error: str,
    error_type: str,
    traceback: str | None = None,
    fatal: bool = False,
) -> Path:
    """Persist table quality failure details alongside existing metadata."""

    path = Path(meta_path)
    metadata = _load_metadata(path)

    failure: dict[str, Any] = {
        "status": "failed",
        "error": error,
        "error_type": error_type,
        "captured_at": datetime.now(UTC).isoformat(),
    }
    if traceback:
        failure["traceback"] = traceback
    if fatal:
        failure["fatal"] = True

    metadata["quality_report"] = failure

    with open_atomic(path, encoding="utf-8") as fh:
        yaml.safe_dump(metadata, fh, sort_keys=False)
        logger.info("metadata_quality_status", path=str(path), status="failed")

    return path
