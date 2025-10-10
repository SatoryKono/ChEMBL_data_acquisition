"""Metadata helpers for generated CSV artefacts."""

from __future__ import annotations

import sys
from collections.abc import Mapping, Sequence
from datetime import UTC as _UTC, datetime as _datetime
from pathlib import Path
from typing import Any

from ..common.metadata_writer import write_meta_yaml as _write_meta_yaml_impl
from ..common.run_context import get_current
from ..config import Config, _serialize_paths

# ---------------------------------------------------------------------------
# Compatibility aliases
# ---------------------------------------------------------------------------

# Historically callers patched :mod:`library.io.metadata.datetime` for
# deterministic timestamps.  Expose ``datetime`` and ``UTC`` to preserve that
# contract for legacy tests and scripts.
UTC = _UTC
datetime = _datetime


def write_meta_yaml(
    path: Path | str,
    cfg: Config | None = None,
    *,
    columns: Sequence[str] | None = None,
    dtypes: Mapping[str, str] | None = None,
    generated_at: str | None = None,
) -> Path:
    """Write metadata for ``path`` to ``<path>.meta.yaml``."""

    if dtypes is None and columns is not None:
        dtypes = {col: "string" for col in columns}

    config_mapping: Mapping[str, Any] | None = None
    if cfg is not None:
        config_mapping = _serialize_paths(cfg.to_dict())

    timestamp = generated_at
    if timestamp is None:
        context = get_current()
        if context is not None and context.generated_at:
            timestamp = context.generated_at
    if timestamp is None:
        timestamp = datetime.now(UTC).isoformat()

    return _write_meta_yaml_impl(
        path,
        command=" ".join(sys.argv),
        config=config_mapping,
        inputs={},
        stats={},
        schema=None,
        generated_at=timestamp,
        columns=columns,
        dtypes=dtypes,
    )


__all__ = ["UTC", "datetime", "write_meta_yaml"]
