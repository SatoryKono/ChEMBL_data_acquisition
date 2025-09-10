"""Load project configuration from YAML files.

This module provides a small helper to load the optional ``config.yaml``
file located at the repository root.  Consumers can supply a custom path if
desired.  Failure to read or parse the configuration results in a
:class:`ConfigError` exception.

The configuration file is expected to be valid YAML.  If the file is missing
or contains malformed content, :func:`load_config` raises :class:`ConfigError`
with a descriptive message.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


class ConfigError(Exception):
    """Raised when the configuration file cannot be loaded."""


# Default location of the configuration file at the project root
DEFAULT_PATH = Path(__file__).resolve().parents[1] / "config.yaml"


def load_config(path: str | Path = DEFAULT_PATH) -> dict[str, Any]:
    """Return configuration settings from ``path``.

    Parameters
    ----------
    path:
        Location of the configuration file. Defaults to the repository's
        ``config.yaml``.

    Returns
    -------
    dict[str, Any]
        Mapping loaded from the YAML document. Empty files result in an empty
        dictionary.

    Raises
    ------
    ConfigError
        If ``path`` does not exist or contains invalid YAML.
    """
    try:
        with Path(path).open("r", encoding="utf8") as fh:
            data = yaml.safe_load(fh) or {}
    except FileNotFoundError as exc:
        raise ConfigError(f"configuration file not found: {path}") from exc
    except yaml.YAMLError as exc:
        raise ConfigError(f"invalid YAML in configuration file: {path}") from exc
    return data
