"""Aggregation utilities for publication metadata."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any

from ..log import logger

TITLE_PREVIEW_LIMIT = 80
TITLE_SUFFIX = "..."

__all__ = ["merge_records", "print_results"]


def merge_records(*records: Mapping[str, Any]) -> dict[str, Any]:
    """Merge mapping objects, later values overriding earlier ones.

    Parameters
    ----------
    *records:
        Dictionaries containing metadata.

    Returns
    -------
    dict
        Combined dictionary containing keys from all inputs.
    """

    merged: dict[str, Any] = {}
    for rec in records:
        merged.update(rec)
    return merged


def _emit_output(log: Any, level: str, output: str) -> None:
    """Send ``output`` to ``log`` using ``level`` when possible."""

    method = getattr(log, level.lower(), None)
    if callable(method):
        try:
            method("pubmed_record_dump", dump=output)
        except TypeError:
            method(output)
        return

    log_method = getattr(log, "log", None)
    if callable(log_method):
        level_name = level.upper()
        try:
            log_method(level_name, "pubmed_record_dump", dump=output)
        except TypeError:
            level_no = getattr(logging, level_name, logging.INFO)
            log_method(level_no, output)


def _first_non_empty(record: Mapping[str, Any], keys: tuple[str, ...]) -> str:
    """Return the first non-empty string value from ``record`` keyed by ``keys``."""

    for key in keys:
        value = record.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def _shorten_title(title: str) -> str:
    """Return a shortened representation of ``title`` suitable for logs."""

    clean = title.strip()
    if len(clean) <= TITLE_PREVIEW_LIMIT:
        return clean
    return clean[: TITLE_PREVIEW_LIMIT - len(TITLE_SUFFIX)] + TITLE_SUFFIX


def print_results(records: list[dict[str, str]], *, level: str = "DEBUG") -> None:
    """Log result records instead of printing to ``stdout``."""

    log = logger
    try:
        from tabulate import tabulate

        use_tabulate = True
    except Exception:
        use_tabulate = False

    display_records = []
    for rec in records:
        pmid = _first_non_empty(
            rec,
            (
                "PubMed.PMID",
                "scholar.PMID",
                "PMID",
            ),
        )
        doi = _first_non_empty(
            rec,
            (
                "PubMed.DOI",
                "scholar.DOI",
                "crossref.DOI",
                "DOI",
            ),
        )
        title = _shorten_title(
            _first_non_empty(
                rec,
                (
                    "PubMed.ArticleTitle",
                    "crossref.Title",
                    "Title",
                ),
            )
        )
        display_records.append({"PMID": pmid, "DOI": doi, "Title": title})

    if use_tabulate:
        output = tabulate(display_records, headers="keys")
    else:
        import json

        output = json.dumps(display_records, ensure_ascii=False, indent=2)

    try:
        _emit_output(log, level, output)
    except UnicodeEncodeError:
        import sys

        encoded = output.encode(sys.stdout.encoding or "utf-8", errors="replace")
        sys.stdout.buffer.write(encoded + b"\n")
