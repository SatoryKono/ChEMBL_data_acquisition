"""Helpers for configuring reusable table-quality reporting hooks."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from functools import partial
from pathlib import Path
from typing import Any

from .table_quality import analyze_table_quality

TableQualitySubject = Any
TableQualityHook = Callable[[TableQualitySubject], None]

__all__ = [
    "TableQualityHook",
    "TableQualitySubject",
    "build_table_quality_hook",
    "is_quality_enabled",
]


def _cfg_value(cfg: Any, key: str, default: Any = None) -> Any:
    if cfg is None:
        return default
    if isinstance(cfg, Mapping):
        return cfg.get(key, default)
    return getattr(cfg, key, default)


def _normalise_columns(columns: Any) -> tuple[str, ...] | None:
    if columns is None:
        return None
    if isinstance(columns, str):
        return (columns,)
    if isinstance(columns, Sequence):
        return tuple(str(column) for column in columns)
    raise TypeError("include/exclude columns must be a sequence of strings or None")


def _normalise_table_name(table_name: str | Path) -> str:
    name = str(table_name)
    if not name:
        raise ValueError("table_name must be a non-empty string")
    return name


def _resolve_destination(
    destination: Path | str | None,
    table_name: str | Path,
) -> Path:
    if destination is not None:
        return Path(destination)
    try:
        path = Path(table_name)
    except (TypeError, ValueError):
        return Path(".")
    parent = path.parent
    return parent if str(parent) else Path(".")


def is_quality_enabled(cfg: Any) -> bool:
    """Return ``True`` if table-quality profiling is enabled for ``cfg``."""

    return bool(_cfg_value(cfg, "enable", False))


def build_table_quality_hook(
    cfg: Any,
    *,
    table_name: str | Path,
    destination: Path | str | None = None,
) -> TableQualityHook:
    """Return a callable that executes :func:`analyze_table_quality` when enabled."""

    if not is_quality_enabled(cfg):
        return _noop_table_quality

    table_name_str = _normalise_table_name(table_name)
    destination_dir = _resolve_destination(destination, table_name)
    include_columns = _normalise_columns(_cfg_value(cfg, "include_columns"))
    exclude_columns = _normalise_columns(_cfg_value(cfg, "exclude_columns"))
    sample_rows = _cfg_value(cfg, "sample_rows")

    return partial(
        analyze_table_quality,
        table_name=table_name_str,
        destination_dir=destination_dir,
        sample_rows=sample_rows,
        include_columns=include_columns,
        exclude_columns=exclude_columns,
    )


def _noop_table_quality(_: TableQualitySubject) -> None:
    return None
