"""Configuration package exposing shared file-system paths."""

from __future__ import annotations

import atexit
from contextlib import ExitStack
from importlib import resources
from pathlib import Path

_STACK = ExitStack()
atexit.register(_STACK.close)
_CONFIG_FILENAME = "config.yaml"


def _resource_path(*parts: str) -> Path:
    traversable = resources.files(__name__)
    for part in parts:
        traversable = traversable.joinpath(part)
    return Path(_STACK.enter_context(resources.as_file(traversable)))


DEFAULT_CONFIG_PATH = _resource_path(_CONFIG_FILENAME)
CONFIG_DIR = DEFAULT_CONFIG_PATH.parent

__all__ = ["CONFIG_DIR", "DEFAULT_CONFIG_PATH"]
