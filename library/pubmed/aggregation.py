"""Aggregation utilities for publication metadata."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any

from ..log import logger

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
        d = rec.copy()
        title = (
            d.get("PubMed.ArticleTitle")
            or d.get("crossref.Title")
            or d.get("Title")
            or ""
        )
        d["Title"] = title[:77] + "..." if len(title) > 80 else title
        display_records.append(d)

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
