"""Metadata helpers for generated CSV artefacts."""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from ..common.run_context import get_current
from ..config import Config, _mask_secrets, _serialize_paths
from ..git_utils import _git_sha
from ..utils.atomic import open_atomic


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

    context = get_current()
    timestamp = generated_at
    if timestamp is None and context is not None and context.generated_at:
        timestamp = context.generated_at
    if timestamp is None:
        timestamp = datetime.now(timezone.utc).isoformat()

    meta = {
        "generated_at": timestamp,
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "columns": list(columns or []),
        "dtypes": dict(dtypes or {}),
        "config": (
            _mask_secrets(_serialize_paths(cfg.to_dict())) if cfg is not None else {}
        ),
    }
    meta_path = Path(f"{path}.meta.yaml")
    with open_atomic(meta_path, encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    meta_lock_path = meta_path.with_name(meta_path.name + ".lock")
    meta_lock_path.unlink(missing_ok=True)
    return meta_path
