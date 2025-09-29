"""Metadata helpers for generated CSV artefacts."""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from ..config import Config, _serialize_paths
from ..git_utils import _git_sha


def write_meta_yaml(
    path: Path | str,
    cfg: Config | None = None,
    *,
    columns: Sequence[str] | None = None,
    dtypes: Mapping[str, str] | None = None,
) -> Path:
    """Write metadata for ``path`` to ``<path>.meta.yaml``."""

    if dtypes is None and columns is not None:
        dtypes = {col: "string" for col in columns}

    meta = {
        "generated_at": datetime.now().isoformat(),
        "git_sha": _git_sha(),
        "command": " ".join(sys.argv),
        "columns": list(columns or []),
        "dtypes": dict(dtypes or {}),
        "config": _serialize_paths(cfg.to_dict()) if cfg is not None else {},
    }
    meta_path = Path(f"{path}.meta.yaml")
    with meta_path.open("w", encoding="utf8") as fh:
        yaml.safe_dump(meta, fh, sort_keys=False)
    return meta_path
