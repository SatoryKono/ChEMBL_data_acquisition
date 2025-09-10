"""Configuration loading utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml


def load_config(path: str | Path) -> Dict[str, Any]:
    """Load a YAML configuration file.

    Parameters
    ----------
    path:
        Location of the configuration file.

    Returns
    -------
    dict[str, Any]
        Parsed configuration settings. An empty dictionary is returned if the
        file does not exist.

    Raises
    ------
    ValueError
        If the file exists but contains invalid YAML.
    """
    cfg_path = Path(path)
    if not cfg_path.exists():
        return {}
    try:
        with cfg_path.open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as exc:  # pragma: no cover - parse errors are rare
        raise ValueError(
            f"failed to parse configuration file {cfg_path}: {exc}"
        ) from exc
    return data or {}
