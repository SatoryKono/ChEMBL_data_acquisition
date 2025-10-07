"""Declarative step callables for assay post-processing."""
from __future__ import annotations

from typing import Mapping, Sequence

from .._utils import build_step_payload

__all__ = ["enrich_metadata", "compute_quality_flags", "export_summary"]


def enrich_metadata(
    *,
    assay_lookup: str,
    fallback_relationship: str,
    include_columns: Sequence[str] | None = None,
) -> Mapping[str, object]:
    """Describe assay metadata enrichment from lookup tables."""

    return build_step_payload(
        "assays.enrich_metadata",
        assay_lookup=assay_lookup,
        fallback_relationship=fallback_relationship,
        include_columns=tuple(include_columns) if include_columns is not None else (),
    )


def compute_quality_flags(
    *,
    quality_ruleset: str,
    strict: bool,
    warn_only: bool,
) -> Mapping[str, object]:
    """Describe computing boolean QA flags for assay exports."""

    return build_step_payload(
        "assays.compute_quality_flags",
        quality_ruleset=quality_ruleset,
        strict=strict,
        warn_only=warn_only,
    )


def export_summary(
    *,
    output_dir: str,
    filename: str,
    include_stats: bool = True,
) -> Mapping[str, object]:
    """Describe exporting aggregate assay statistics."""

    return build_step_payload(
        "assays.export_summary",
        output_dir=output_dir,
        filename=filename,
        include_stats=include_stats,
    )
