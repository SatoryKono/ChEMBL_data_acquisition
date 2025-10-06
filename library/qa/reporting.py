"""Helpers for configuring table-quality reporting.

This module centralises construction of table-quality hooks used across
pipelines.  Individual scripts configure a :class:`TableQualityReporter`
instance and then either build a hook compatible with
``library.cli_utils.run_pipeline`` or run the analysis directly on cached
profilers.  Consolidating this logic ensures all entry points share the
same configuration defaults (for example sampling options and the
``fatal_on_error`` flag).
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from .table_quality import TableQualityProfiler, analyze_table_quality

if TYPE_CHECKING:  # pragma: no cover - imported for type checking only
    from library.config import Config, DocQualityCfg


TableQualityHook = Callable[[Path], None]
TableQualitySource = (
    pd.DataFrame | str | Path | Iterable[pd.DataFrame] | TableQualityProfiler
)


def _as_doc_quality_cfg(
    cfg: Config | DocQualityCfg | None,
) -> DocQualityCfg | None:
    """Return the ``DocQualityCfg`` instance extracted from ``cfg``."""

    if cfg is None:
        return None

    if hasattr(cfg, "enable") and hasattr(cfg, "sample_rows"):
        # ``cfg`` already behaves like ``DocQualityCfg``
        return cfg  # type: ignore[return-value]

    system_cfg = getattr(cfg, "system", None)
    if system_cfg is None:
        return None
    return getattr(system_cfg, "doc_quality", None)


def _normalise_columns(columns: Iterable[str] | None) -> tuple[str, ...] | None:
    if columns is None:
        return None
    return tuple(columns)


@dataclass(frozen=True)
class TableQualityReporter:
    """Configure and execute table-quality analysis."""

    enabled: bool = False
    sample_rows: int | None = None
    include_columns: tuple[str, ...] | None = None
    exclude_columns: tuple[str, ...] | None = None
    fatal_on_error: bool = False

    @classmethod
    def from_config(
        cls, cfg: Config | DocQualityCfg | None
    ) -> "TableQualityReporter":
        doc_quality_cfg = _as_doc_quality_cfg(cfg)
        if doc_quality_cfg is None:
            return cls()
        return cls(
            enabled=bool(getattr(doc_quality_cfg, "enable", False)),
            sample_rows=getattr(doc_quality_cfg, "sample_rows", None),
            include_columns=_normalise_columns(
                getattr(doc_quality_cfg, "include_columns", None)
            ),
            exclude_columns=_normalise_columns(
                getattr(doc_quality_cfg, "exclude_columns", None)
            ),
            fatal_on_error=bool(getattr(doc_quality_cfg, "fatal_on_error", False)),
        )

    def build_hook(
        self, *, table_name: str, destination_dir: Path | str
    ) -> TableQualityHook:
        """Return a callable suitable for ``run_pipeline`` hooks."""

        if not self.enabled:
            return _noop_quality_hook

        resolved_name = str(table_name)
        destination = Path(destination_dir)

        def _hook(path: Path) -> None:
            analyze_table_quality(
                path,
                table_name=resolved_name,
                destination_dir=destination,
                sample_rows=self.sample_rows,
                include_columns=self.include_columns,
                exclude_columns=self.exclude_columns,
            )

        return _hook

    def run(
        self,
        table: TableQualitySource,
        *,
        table_name: str,
        destination_dir: Path | str,
    ) -> tuple[pd.DataFrame, pd.DataFrame] | None:
        """Execute table-quality analysis when enabled."""

        if not self.enabled:
            return None

        return analyze_table_quality(
            table,
            table_name=table_name,
            destination_dir=destination_dir,
            sample_rows=self.sample_rows,
            include_columns=self.include_columns,
            exclude_columns=self.exclude_columns,
        )


def _noop_quality_hook(_: Path) -> None:
    return None


__all__ = ["TableQualityHook", "TableQualityReporter", "TableQualitySource"]
