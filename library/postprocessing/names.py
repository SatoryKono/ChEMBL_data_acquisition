"""Helpers for generating target name exports.

The helper mirrors the legacy Power Query workbook that produced auxiliary
name tables for reporting.  It extracts textual identifiers from the merged
targets export, normalises them and emits a long-form table where each row
represents a distinct name attributed to a target.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from . import helpers
from ..common.log import logger

# Columns in the targets table used to derive names.  Each entry maps the column
# name to a descriptive label stored alongside the extracted token so consumers
# can identify the provenance of a particular name.
NAME_COLUMN_SOURCES: tuple[tuple[str, str], ...] = (
    ("pref_name", "chembl_preferred"),
    ("protein_name_canonical", "uniprot_canonical"),
    ("protein_name_alt", "uniprot_alternative"),
    ("synonyms", "chembl_synonym"),
    ("protein_synonym_list", "uniprot_synonym"),
    ("gtop_synonyms", "gtop_synonym"),
    ("gene_symbol", "gene_symbol"),
    ("gene_symbol_list", "gene_symbol"),
    ("isoform_names", "isoform_name"),
    ("isoform_synonyms", "isoform_synonym"),
    ("recommendedName", "uniprot_recommended"),
    ("secondaryAccessionNames", "uniprot_secondary"),
)

# Columns that should be split on the pipe delimiter.  Other columns are treated
# as single textual values even if a pipe character is present.
PIPE_SPLIT_COLUMNS: frozenset[str] = frozenset(
    {
        "synonyms",
        "protein_synonym_list",
        "gtop_synonyms",
        "gene_symbol_list",
        "isoform_names",
        "isoform_synonyms",
        "secondaryAccessionNames",
    }
)

# Canonical column order for the emitted names table.
TARGET_NAMES_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "uniprot_id_primary",
    "name",
    "name_type",
    "source_column",
)

# Sentinel values considered as missing tokens when splitting pipe-delimited
# strings.  The list mirrors guards in the legacy workbook.
EMPTY_TOKEN_MARKERS: frozenset[str] = frozenset({"", "-", "n/a", "none"})


def _iter_name_tokens(value: object, *, split: bool) -> Iterable[str]:
    """Yield normalised name tokens extracted from ``value``."""

    text = helpers.normalise_text(value)
    if not text:
        return []

    if not split:
        return [text]

    tokens: list[str] = []
    for raw in text.split("|"):
        token = helpers.normalise_text(raw)
        if not token:
            continue
        lower = token.lower()
        if lower in EMPTY_TOKEN_MARKERS:
            continue
        tokens.append(token)
    return tokens


def _build_names_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Project ``frame`` onto the long-form target names table."""

    records: list[dict[str, Any]] = []
    column_set = set(frame.columns)

    for idx, row in frame.iterrows():
        target_id = helpers.normalise_text(row.get("target_chembl_id"))
        if not target_id:
            continue
        uniprot_id = helpers.normalise_text(row.get("uniprot_id_primary"))
        for column, label in NAME_COLUMN_SOURCES:
            if column not in column_set:
                continue
            split = column in PIPE_SPLIT_COLUMNS
            for token in _iter_name_tokens(row.get(column), split=split):
                records.append(
                    {
                        "target_chembl_id": target_id,
                        "uniprot_id_primary": uniprot_id,
                        "name": token,
                        "name_type": label,
                        "source_column": column,
                    }
                )

    names_df = pd.DataFrame.from_records(records, columns=TARGET_NAMES_COLUMNS)
    if names_df.empty:
        return pd.DataFrame(columns=TARGET_NAMES_COLUMNS, dtype="string")

    names_df = helpers.ensure_string_columns(names_df, TARGET_NAMES_COLUMNS)
    names_df = names_df.drop_duplicates().sort_values(
        by=["target_chembl_id", "name", "name_type"], kind="mergesort"
    )
    names_df = names_df.reset_index(drop=True)
    return names_df


def _summarise_contrion(series: pd.Series | None) -> dict[str, int]:
    """Return counts for the ``contrion`` column if present."""

    if series is None:
        return {"contrion_unique": 0, "contrion_non_empty": 0, "contrion_total": 0}

    tokens: set[str] = set()
    non_empty = 0
    total = 0
    for value in series.dropna().astype(str):
        text = helpers.normalise_text(value)
        if not text or text in EMPTY_TOKEN_MARKERS:
            continue
        parts = [helpers.normalise_text(part) for part in text.split("|")]
        parts = [part for part in parts if part and part not in EMPTY_TOKEN_MARKERS]
        if not parts:
            continue
        non_empty += 1
        total += len(parts)
        tokens.update(parts)
    return {
        "contrion_unique": len(tokens),
        "contrion_non_empty": non_empty,
        "contrion_total": total,
    }


def _summarise_active_component_type(series: pd.Series | None) -> dict[str, int]:
    """Return counts per ``active_component_type`` value."""

    if series is None:
        return {}

    normalised = (
        series.fillna("<NA>")
        .astype(str)
        .map(lambda value: helpers.normalise_text(value) or "<EMPTY>")
    )
    counts = normalised.value_counts(dropna=False)
    return {str(key): int(value) for key, value in counts.items()}


def process_target_names(input_path: str | Path, *, verbose: bool = False) -> dict[str, Any]:
    """Process ``input_path`` and emit the target names table."""

    source_path = Path(input_path)
    frame = helpers.read_csv_with_fallbacks(source_path)
    frame = helpers.ensure_string_columns(frame, frame.columns)

    names_df = _build_names_table(frame)
    output_path = source_path.with_name(f"names.{source_path.name}")
    helpers.write_csv(names_df, output_path, columns=TARGET_NAMES_COLUMNS)

    summary: dict[str, Any] = {
        "rows_before": int(len(frame)),
        "rows_after": int(len(names_df)),
    }
    contrion_summary = _summarise_contrion(frame.get("contrion"))
    summary.update(contrion_summary)
    active_summary = _summarise_active_component_type(frame.get("active_component_type"))
    if active_summary:
        summary["active_component_type"] = active_summary
    if verbose:
        logger.info(
            "target_names_helper_done",
            path=str(output_path),
            rows=len(names_df),
        )
    return {"path": str(output_path), "summary": summary}
