"""Shared helpers for writing metadata sidecar files.

This module centralises the implementation previously duplicated between the
``library.common.metadata`` and ``library.io.metadata`` modules.  The
``write_meta_yaml`` entry point supports the full set of rich metadata fields
used by the pipelines while remaining backwards compatible with the simpler
``library.io.metadata`` shim.
"""

from __future__ import annotations

import hashlib
import platform
import shlex
import sys
from collections.abc import Mapping, Sequence
from datetime import datetime, timedelta, timezone
from functools import lru_cache
from pathlib import Path
from typing import Any, TypedDict

import yaml

from config.paths import DICTIONARY_DIR

from ..config import _mask_secrets
from ..project_version import get_pipeline_version
from ..resources.dictionaries import DictionaryManifestError, get_resource
from ..utils.atomic import open_atomic
from .git import _git_sha
from .log import logger
from .run_context import get_current

# ``datetime.UTC`` is only available in Python 3.11 and later.  Using
# ``timezone.utc`` keeps the implementation compatible with older runtimes.
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
    pubchem_augmentation_enabled: bool
    pubchem_fallback_applied: bool
    missing_ids_sample: list[str]
    missing_molecule_ids: list[str]
    missing_molecule_ids_total: int
    missing_molecule_ids_truncated: bool
    missing_molecule_ids_count: int
    chunk_fetch_failure_chunks: int
    chunk_fetch_failure_ids: list[str]
    chunk_fetch_failure_ids_total: int
    chunk_fetch_failure_ids_truncated: bool


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


def file_sha256(path: Path | str) -> str:
    """Return the SHA256 checksum of ``path``."""

    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_DETERMINISTIC_TIMESTAMP_BASE = datetime(2000, 1, 1, tzinfo=UTC)
_DETERMINISTIC_TIMESTAMP_WINDOW_SECONDS = 100 * 365 * 24 * 60 * 60


def compute_generated_at(
    path: Path | str,
    *,
    command: str,
    invocation: Sequence[str] | None = None,
) -> str:
    """Return a deterministic timestamp derived from the invocation context."""

    key_parts = [str(Path(path)), command]
    if invocation:
        key_parts.extend(invocation)
    key = "\u241f".join(key_parts).encode("utf-8")

    digest = hashlib.sha256(key).digest()
    seconds = int.from_bytes(digest[:6], "big") % _DETERMINISTIC_TIMESTAMP_WINDOW_SECONDS
    microseconds = int.from_bytes(digest[6:9], "big") % 1_000_000
    timestamp = _DETERMINISTIC_TIMESTAMP_BASE + timedelta(
        seconds=seconds,
        microseconds=microseconds,
    )
    return timestamp.isoformat()


def write_meta_yaml(
    csv_path: Path | str,
    *,
    command: str | None = None,
    config: Mapping[str, Any] | None = None,
    inputs: Mapping[str, Any] | None = None,
    stats: Mapping[str, Any] | None = None,
    schema: str | None = None,
    invocation: Sequence[str] | None = None,
    extra_metadata: Mapping[str, Any] | None = None,
    dictionary_resources: Sequence[str] | None = None,
    generated_at: str | None = None,
    allow_nondeterministic_timestamp: bool = False,
    columns: Sequence[str] | None = None,
    dtypes: Mapping[str, str] | None = None,
) -> Path:
    """Write metadata for ``csv_path`` to ``<csv_path>.meta.yaml``."""

    path = Path(csv_path)
    meta_path = path.with_name(path.name + ".meta.yaml")

    existing = _load_metadata(meta_path)

    metadata: dict[str, Any] = dict(existing)
    context = get_current()

    command_tokens = _split_command(command)
    if command_tokens is not None:
        normalised_tokens = _normalise_command_tokens(command_tokens, output_path=path)
        command_str = " ".join(normalised_tokens)
    else:
        command_str = command if command is not None else " ".join(sys.argv)
        normalised_tokens = None

    invocation_for_seed: Sequence[str] | None
    if invocation is not None:
        invocation_for_seed = list(invocation)
    else:
        existing_invocation = existing.get("invocation")
        if isinstance(existing_invocation, Sequence) and not isinstance(
            existing_invocation,
            (str, bytes),
        ):
            invocation_for_seed = list(existing_invocation)
        else:
            invocation_for_seed = None

    timestamp = generated_at
    if timestamp is None and context is not None and context.generated_at:
        timestamp = context.generated_at
    if timestamp is None:
        if allow_nondeterministic_timestamp:
            timestamp = datetime.now(UTC).isoformat()
        else:
            timestamp = compute_generated_at(
                path,
                command=command_str,
                invocation=invocation_for_seed,
            )

    metadata.update(
        {
            "generated_at": timestamp,
            "git_sha": _git_sha(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "command": command_str,
            "config": _mask_secrets(dict(config or {})),
            "inputs": dict(inputs or {}),
            "stats": dict(stats or {}),
            "schema": schema,
            "columns": list(columns or []),
            "dtypes": dict(dtypes or {}),
            "pipeline_version": get_pipeline_version(),
            "output_path": str(path),
        }
    )

    if normalised_tokens is not None:
        metadata["command_args"] = list(normalised_tokens)
    elif "command_args" in metadata:
        metadata.pop("command_args", None)

    if invocation is not None:
        metadata["invocation"] = list(invocation)
    elif "invocation" in metadata and not metadata["invocation"]:
        # Normalise empty invocation lists that may be persisted in pre-existing
        # metadata files.
        metadata.pop("invocation", None)

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


_COMMAND_OUTPUT_PLACEHOLDER = "<OUTPUT_PATH>"
_COMMAND_ABS_PATH_PLACEHOLDER = "<ABS_PATH>"


def _split_command(command: str | None) -> Sequence[str] | None:
    """Return tokens extracted from ``command`` when possible."""

    if command is None:
        return list(sys.argv)

    try:
        return shlex.split(command)
    except ValueError:
        logger.debug("metadata_writer_command_split_failed", command=command)
        return None


def _normalise_command_tokens(
    tokens: Sequence[str],
    *,
    output_path: Path,
) -> list[str]:
    """Return ``tokens`` with absolute paths replaced by placeholders."""

    normalised: list[str] = []
    skip_next = False
    output_str = str(output_path)

    for index, token in enumerate(tokens):
        if skip_next:
            skip_next = False
            continue

        if token == "--final-out":
            normalised.append(token)
            if index + 1 < len(tokens):
                normalised.append(_COMMAND_OUTPUT_PLACEHOLDER)
                skip_next = True
            continue

        if token.startswith("--final-out="):
            prefix, _, _ = token.partition("=")
            normalised.append(f"{prefix}={_COMMAND_OUTPUT_PLACEHOLDER}")
            continue

        if token == output_str:
            normalised.append(_COMMAND_OUTPUT_PLACEHOLDER)
            continue

        name, sep, value = token.partition("=")
        if sep and value.startswith("/"):
            normalised.append(f"{name}={_COMMAND_ABS_PATH_PLACEHOLDER}")
            continue

        if token.startswith("/"):
            normalised.append(_COMMAND_ABS_PATH_PLACEHOLDER)
            continue

        normalised.append(token)

    return normalised


__all__ = [
    "Stats",
    "file_sha256",
    "record_quality_failure",
    "compute_generated_at",
    "write_meta_yaml",
]

