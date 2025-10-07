"""Declarative step callables for document post-processing."""
from __future__ import annotations

from typing import Mapping, Sequence

from .._utils import build_step_payload

__all__ = ["normalize_documents", "classify_types", "export_quality_report"]


def normalize_documents(
    *,
    input_csv: str,
    text_columns: Sequence[str],
    strip_html: bool,
) -> Mapping[str, object]:
    """Describe normalising textual document columns."""

    return build_step_payload(
        "documents.normalize_documents",
        input_csv=input_csv,
        text_columns=tuple(text_columns),
        strip_html=strip_html,
    )


def classify_types(
    *,
    classifier_model: str,
    review_threshold: float,
    experimental_threshold: float,
) -> Mapping[str, object]:
    """Describe applying document type classification thresholds."""

    return build_step_payload(
        "documents.classify_types",
        classifier_model=classifier_model,
        review_threshold=review_threshold,
        experimental_threshold=experimental_threshold,
    )


def export_quality_report(
    *,
    output_dir: str,
    filename: str,
    include_logs: bool,
) -> Mapping[str, object]:
    """Describe exporting QA reports for document tables."""

    return build_step_payload(
        "documents.export_quality_report",
        output_dir=output_dir,
        filename=filename,
        include_logs=include_logs,
    )
