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
from functools import lru_cache
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, TypedDict

import yaml

from config.paths import DICTIONARY_DIR

from ..config import _mask_secrets
from ..resources.dictionaries import DictionaryManifestError, get_resource
from .git import _git_sha
from .log import logger
from ..utils.atomic import open_atomic
from ..project_version import get_pipeline_version


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
    missing_molecule_ids: list[str]
    missing_molecule_ids_count: int
    chunk_fetch_failure_chunks: int
    chunk_fetch_failure_ids: list[str]
    chunk_fetch_failure_ids_total: int
    chunk_fetch_failure_ids_truncated: bool


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
    dictionary_resources: Sequence[str] | None = None,
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
    dictionary_resources:
        Optional sequence of dictionary resource identifiers declared in the
        bundled manifest.  When provided, the metadata output is enriched with
        the corresponding version and SHA256 checksum for every listed
        resource.

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
            "pipeline_version": get_pipeline_version(),
        }
    )
    if invocation:
        metadata["invocation"] = list(invocation)
    if extra_metadata:
        metadata.update(dict(extra_metadata))

    if dictionary_resources:
        manifest_metadata = _dictionary_manifest_metadata()
        dictionaries: dict[str, Mapping[str, str]] = {}
        for name in dict.fromkeys(dictionary_resources):
            try:
                resource = get_resource(name)
            except KeyError as exc:
                logger.warning(
                    "metadata_dictionary_lookup_failed",
                    resource=name,
                    error=str(exc),
                )
                entry = manifest_metadata.get(name)
                if entry:
                    dictionaries[name] = dict(entry)
                continue
            except DictionaryManifestError as exc:
                logger.warning(
                    "metadata_dictionary_lookup_failed",
                    resource=name,
                    error=str(exc),
                )
                entry = manifest_metadata.get(name)
                if entry:
                    dictionaries[name] = dict(entry)
                continue
            dictionaries[name] = {
                "version": resource.version,
                "sha256": resource.sha256,
            }
        if dictionaries:
            existing_dictionaries = metadata.get("dictionaries")
            if isinstance(existing_dictionaries, Mapping):
                merged = dict(existing_dictionaries)
                merged.update(dictionaries)
                metadata["dictionaries"] = merged
            else:
                metadata["dictionaries"] = dictionaries

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
@lru_cache(maxsize=1)
def _dictionary_manifest_metadata() -> Mapping[str, Mapping[str, str]]:
    manifest_path = DICTIONARY_DIR / "manifest.yaml"
    try:
        with manifest_path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    except OSError:
        return {}

    resources = data.get("resources")
    if not isinstance(resources, Mapping):
        return {}

    entries: dict[str, Mapping[str, str]] = {}
    for name, meta in resources.items():
        if not isinstance(meta, Mapping):
            continue
        version = meta.get("version")
        sha256 = meta.get("sha256")
        if isinstance(version, str) and isinstance(sha256, str):
            entries[name] = {"version": version, "sha256": sha256}
    return entries

