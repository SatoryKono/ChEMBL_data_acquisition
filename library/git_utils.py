"""Compatibility layer for Git helpers.

Historically :mod:`library.git_utils` contained its own copy of the Git helper
logic.  Maintaining two separate implementations led to duplicated behaviour
and—importantly for Windows users—duplicated side effects such as logging
``git_sha_fallback`` twice for a single pipeline run.  The canonical
implementation now lives in :mod:`library.common.git`; this module simply
re-exports the public helpers so existing imports keep working while sharing a
single cached implementation.
"""

from __future__ import annotations

from .common.git import (
    _ensure_text,
    _error_payload,
    _format_error,
    _format_subprocess_error,
    _git_sha,
    _is_sha,
    _log_fallback,
    _normalise_returncode,
    _read_head_sha,
    _read_packed_ref,
    _repo_root,
    _resolve_git_dir,
    _serialise_cmd,
)

__all__ = [
    "_ensure_text",
    "_error_payload",
    "_format_error",
    "_format_subprocess_error",
    "_git_sha",
    "_is_sha",
    "_log_fallback",
    "_normalise_returncode",
    "_read_head_sha",
    "_read_packed_ref",
    "_repo_root",
    "_resolve_git_dir",
    "_serialise_cmd",
]
