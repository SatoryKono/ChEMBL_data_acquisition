"""Declarative step callables for the activity post-processing pipeline."""
from __future__ import annotations

from typing import Mapping

from .._utils import build_step_payload

__all__ = ["compute_bounds", "annotate_targets", "export_extended"]


def compute_bounds(*, input_csv: str, rounding_digits: int = 3, clamp_nonnegative: bool = True) -> Mapping[str, object]:
    """Describe the activity bounds computation step."""

    return build_step_payload(
        "activities.compute_bounds",
        input_csv=input_csv,
        rounding_digits=rounding_digits,
        clamp_nonnegative=clamp_nonnegative,
    )


def annotate_targets(
    *,
    dictionary_csv: str,
    fallback_label: str,
    minimum_confidence: float,
) -> Mapping[str, object]:
    """Describe enrichment of activity rows with target annotations."""

    return build_step_payload(
        "activities.annotate_targets",
        dictionary_csv=dictionary_csv,
        fallback_label=fallback_label,
        minimum_confidence=minimum_confidence,
    )


def export_extended(
    *,
    search_dir: str,
    output_dir: str,
    dictionary_dir: str,
    timestamp_column: str = "timestamp_utc",
) -> Mapping[str, object]:
    """Describe exporting the extended activity artefact."""

    return build_step_payload(
        "activities.export_extended",
        search_dir=search_dir,
        output_dir=output_dir,
        dictionary_dir=dictionary_dir,
        timestamp_column=timestamp_column,
    )
