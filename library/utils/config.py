"""Compatibility helpers for legacy configuration imports.

Historically the command line tooling exposed configuration utilities from
``library.utils.config``.  The configuration package was later reorganised
under :mod:`library.config`, but some scripts and third-party integrations still
import the old module path.  This shim re-exports the public API from the new
location so that existing code continues to work.
"""

from __future__ import annotations

from ..config import (
    Config,
    ConfigError,
    ConfigLoaderError,
    ConfigMetadata,
    ConfigSource,
    ensure_dirs,
    load_config,
    load_yaml_config,
    print_config,
    resolve_config_path,
)
from ..config.loader import DEFAULT_CONFIG_PATH

__all__ = [
    "Config",
    "ConfigError",
    "ConfigLoaderError",
    "ConfigMetadata",
    "ConfigSource",
    "DEFAULT_CONFIG_PATH",
    "ensure_dirs",
    "load_config",
    "load_yaml_config",
    "print_config",
    "resolve_config_path",
]
