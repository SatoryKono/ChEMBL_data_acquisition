"""Declarative step callables for target post-processing."""
from __future__ import annotations

from typing import Mapping, Sequence

from .._utils import build_step_payload

__all__ = ["merge_taxonomy", "prepare_isoforms", "emit_helper_exports"]


def merge_taxonomy(
    *,
    taxonomy_source: str,
    lineage_columns: Sequence[str],
) -> Mapping[str, object]:
    """Describe aligning target taxonomy helpers."""

    return build_step_payload(
        "targets.merge_taxonomy",
        taxonomy_source=taxonomy_source,
        lineage_columns=tuple(lineage_columns),
    )


def prepare_isoforms(
    *,
    isoform_export: str,
    dedupe: bool,
    sort_key: str,
) -> Mapping[str, object]:
    """Describe preparing isoform projections."""

    return build_step_payload(
        "targets.prepare_isoforms",
        isoform_export=isoform_export,
        dedupe=dedupe,
        sort_key=sort_key,
    )


def emit_helper_exports(
    *,
    output_dir: str,
    include_files: Sequence[str],
    timestamp_column: str,
) -> Mapping[str, object]:
    """Describe emitting auxiliary helper exports for targets."""

    return build_step_payload(
        "targets.emit_helper_exports",
        output_dir=output_dir,
        include_files=tuple(include_files),
        timestamp_column=timestamp_column,
    )
