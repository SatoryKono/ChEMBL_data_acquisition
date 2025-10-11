"""Runtime helpers for sharing postprocessing configuration."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator


@dataclass
class _RuntimeState:
    """Container storing runtime configuration shared across modules."""

    export_root: Path | None = None


_STATE = _RuntimeState()


def set_default_export_root(path: Path | str | None) -> None:
    """Set the default directory containing exported CSV artifacts."""

    if path is None:
        _STATE.export_root = None
    else:
        _STATE.export_root = Path(path)


def get_default_export_root() -> Path | None:
    """Return the configured export root directory if available."""

    return _STATE.export_root


@contextmanager
def override_default_export_root(path: Path | str | None) -> Iterator[None]:
    """Temporarily override the configured export root directory."""

    previous = _STATE.export_root
    set_default_export_root(path)
    try:
        yield
    finally:
        _STATE.export_root = previous


def configure_runtime_paths_from_config(cfg: Any) -> None:
    """Populate runtime defaults from a loaded configuration object."""

    io_cfg = getattr(cfg, "io", None)
    if io_cfg is None:
        return

    output_dir = getattr(io_cfg, "output_dir", None)
    if output_dir is None:
        return

    set_default_export_root(output_dir)


__all__ = [
    "configure_runtime_paths_from_config",
    "get_default_export_root",
    "override_default_export_root",
    "set_default_export_root",
]

