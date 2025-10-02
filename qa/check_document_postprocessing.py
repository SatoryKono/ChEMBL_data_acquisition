"""Quality assurance for document post-processing outputs.

Changelog
~~~~~~~~
- 2025-02-??: Initial implementation mirroring legacy Power Query QA.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd

from library.csv_utils import write_csv_deterministic
from library.postprocessing.document import to_text


# ===== Parameters ============================================================
CP1252_ENCODING = "cp1252"
UTF8_ENCODING = "utf-8"
CSV_DELIMITER = ","
DEFAULT_MAX_DIFF_ROWS = 100
KEY_COLUMNS: tuple[str, str, str] = ("PMID", "document_chembl_id", "completed")
REPORT_PREFIX = "qa_document_postprocessing_report_"
DIFF_PREFIX = "qa_document_postprocessing_diff_"
REPORT_JSON_SUFFIX = ".json"
REPORT_MD_SUFFIX = ".md"
DIFF_SUFFIX = ".csv"
DATE_PATTERN = re.compile(r"(20\d{6})")
SUCCESS_STATUS = "passed"
FAILURE_STATUS = "failed"


# ===== Data classes ==========================================================
@dataclass
class QaResult:
    """Structured result emitted by :func:`run_document_postprocessing_check`."""

    passed: bool
    date_code: str
    metrics: Mapping[str, Any]
    report_json: Path
    report_markdown: Path
    diff_csv: Path | None


# ===== Helpers ===============================================================
def _resolve_relative(base_path: Path, relative: str | os.PathLike[str]) -> Path:
    """Return ``relative`` resolved from ``base_path`` honouring Windows paths."""

    relative_text = str(relative)
    candidate = Path(relative_text.replace("\\", os.sep))
    if candidate.is_absolute():
        return candidate
    return (base_path / candidate).resolve()


def _ensure_columns(frame: pd.DataFrame, columns: Iterable[str]) -> pd.DataFrame:
    """Return ``frame`` with ``columns`` ensured (missing ones filled with empty)."""

    result = frame.copy()
    for column in columns:
        if column not in result.columns:
            result[column] = ""
    return result


def _canonicalise(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` with textual normalisation mirroring Power Query."""

    canonical = pd.DataFrame(index=frame.index)
    for column in frame.columns:
        canonical[column] = frame[column].map(to_text)
    return canonical


def _build_index(frame: pd.DataFrame, key_columns: Sequence[str]) -> pd.MultiIndex:
    """Return a MultiIndex composed from ``key_columns`` after canonicalisation."""

    keys = [frame[column].map(to_text) for column in key_columns]
    tuples = list(zip(*keys, strict=False))
    if not tuples:
        return pd.MultiIndex.from_tuples([], names=list(key_columns))
    return pd.MultiIndex.from_tuples(tuples, names=list(key_columns))


def _extract_date_code(path: Path, explicit: str | None = None) -> str:
    """Return the date code associated with ``path`` or ``explicit`` fallback."""

    if explicit:
        return explicit
    match = DATE_PATTERN.search(path.name)
    if match:
        return match.group(1)
    return datetime.utcnow().strftime("%Y%m%d")


def _summarise_missing(index_a: pd.Index, index_b: pd.Index) -> dict[str, Any]:
    """Return counts of elements present in ``index_a`` but not ``index_b``."""

    missing = index_a.difference(index_b)
    sample = [tuple(item) if isinstance(item, tuple) else item for item in missing[:5]]
    return {"count": int(len(missing)), "sample": sample}


def _differences_to_frame(
    mask: pd.DataFrame,
    expected: pd.DataFrame,
    actual: pd.DataFrame,
    *,
    key_columns: Sequence[str],
    limit: int,
) -> pd.DataFrame:
    """Return a tidy DataFrame describing mismatched cells limited by ``limit``."""

    records: list[dict[str, Any]] = []
    columns = list(mask.columns)
    for row_idx, column_idx in zip(*mask.to_numpy().nonzero(), strict=False):
        index_value = mask.index[row_idx]
        if isinstance(index_value, tuple):
            key_map = dict(zip(key_columns, index_value, strict=False))
        else:
            key_map = {key_columns[0]: index_value}
            for extra in key_columns[1:]:
                key_map.setdefault(extra, "")
        record = {
            **key_map,
            "column": columns[column_idx],
            "expected": to_text(expected.iloc[row_idx, column_idx]),
            "actual": to_text(actual.iloc[row_idx, column_idx]),
        }
        records.append(record)
        if len(records) >= limit:
            break
    return pd.DataFrame.from_records(records, columns=[*key_columns, "column", "expected", "actual"])


def _differences_by_column(mask: pd.DataFrame) -> dict[str, int]:
    """Return a dictionary counting mismatches per column."""

    return {column: int(mask[column].sum()) for column in mask.columns}


def _render_markdown(metrics: Mapping[str, Any], diff_path: Path | None) -> str:
    """Return a Markdown summary of the QA run."""

    lines = ["# Document post-processing QA report", ""]
    lines.append(f"*Status*: **{metrics['status']}**")
    lines.append(f"*Reference rows*: {metrics['reference']['rows']}")
    lines.append(f"*Candidate rows*: {metrics['candidate']['rows']}")
    lines.append(
        f"*Cells different*: {metrics['differences']['cells_different']} / "
        f"{metrics['differences']['cells_total']}"
    )
    lines.append(f"*Rows with differences*: {metrics['differences']['rows_with_differences']}")
    lines.append(f"*Issues detected*: {len(metrics['issues'])}")
    if metrics["issues"]:
        lines.append("")
        lines.append("## Issues")
        for issue in metrics["issues"]:
            lines.append(f"- {issue}")
    lines.append("")
    lines.append("## Differences by column")
    lines.append("| Column | Differences |")
    lines.append("| --- | ---: |")
    for column, count in metrics["differences"]["by_column"].items():
        lines.append(f"| {column} | {count} |")
    if diff_path is not None:
        lines.append("")
        lines.append(f"Diff excerpt written to `{diff_path.name}`.")
    return "\n".join(lines) + "\n"


# ===== Core functionality ====================================================
def run_document_postprocessing_check(
    *,
    base_path: str | os.PathLike[str],
    reference_path: str | os.PathLike[str],
    candidate_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str] | None = None,
    date_code: str | None = None,
    delimiter: str = CSV_DELIMITER,
    reference_encoding: str = CP1252_ENCODING,
    candidate_encoding: str = UTF8_ENCODING,
    max_diff_rows: int = DEFAULT_MAX_DIFF_ROWS,
    key_columns: Sequence[str] = KEY_COLUMNS,
) -> QaResult:
    """Compare Power Query reference output with the Python post-processing result."""

    if max_diff_rows <= 0:
        msg = "max_diff_rows must be positive"
        raise ValueError(msg)

    base_dir = Path(base_path).resolve()
    reference_resolved = _resolve_relative(base_dir, reference_path)
    candidate_resolved = _resolve_relative(base_dir, candidate_path)

    if output_dir is None:
        output_directory = candidate_resolved.parent
    else:
        output_directory = _resolve_relative(base_dir, output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    reference_df = pd.read_csv(
        reference_resolved,
        dtype="string",
        encoding=reference_encoding,
        sep=delimiter,
        keep_default_na=False,
    )
    candidate_df = pd.read_csv(
        candidate_resolved,
        dtype="string",
        encoding=candidate_encoding,
        sep=delimiter,
        keep_default_na=False,
    )

    reference_df = _ensure_columns(reference_df, key_columns)
    candidate_df = _ensure_columns(candidate_df, key_columns)

    reference_canonical = _canonicalise(reference_df)
    candidate_canonical = _canonicalise(candidate_df)

    reference_index = _build_index(reference_canonical, key_columns)
    candidate_index = _build_index(candidate_canonical, key_columns)

    reference_canonical.index = reference_index
    candidate_canonical.index = candidate_index

    all_columns = sorted(set(reference_canonical.columns) | set(candidate_canonical.columns))
    reference_canonical = reference_canonical.reindex(columns=all_columns).fillna("")
    candidate_canonical = candidate_canonical.reindex(columns=all_columns).fillna("")

    combined_index = reference_index.union(candidate_index)
    reference_aligned = reference_canonical.reindex(combined_index, fill_value="")
    candidate_aligned = candidate_canonical.reindex(combined_index, fill_value="")

    mask = reference_aligned.ne(candidate_aligned)
    cells_different = int(mask.to_numpy().sum())
    total_cells = int(mask.size)
    rows_with_differences = int(mask.any(axis=1).sum())
    by_column = _differences_by_column(mask)

    missing_in_candidate = _summarise_missing(reference_index, candidate_index)
    missing_in_reference = _summarise_missing(candidate_index, reference_index)

    reference_duplicates = int(reference_index.duplicated().sum())
    candidate_duplicates = int(candidate_index.duplicated().sum())

    missing_reference_columns = sorted(set(all_columns) - set(candidate_df.columns))
    missing_candidate_columns = sorted(set(all_columns) - set(reference_df.columns))

    issues: list[str] = []
    if cells_different:
        issues.append(f"Detected {cells_different} mismatched cells")
    if missing_in_candidate["count"]:
        issues.append(
            f"{missing_in_candidate['count']} keys missing in candidate output"
        )
    if missing_in_reference["count"]:
        issues.append(
            f"{missing_in_reference['count']} keys missing in reference output"
        )
    if reference_duplicates:
        issues.append(f"Reference contains {reference_duplicates} duplicate key rows")
    if candidate_duplicates:
        issues.append(f"Candidate contains {candidate_duplicates} duplicate key rows")
    if missing_reference_columns:
        issues.append(
            "Reference missing columns present in candidate: "
            + ", ".join(missing_reference_columns)
        )
    if missing_candidate_columns:
        issues.append(
            "Candidate missing columns present in reference: "
            + ", ".join(missing_candidate_columns)
        )

    status = SUCCESS_STATUS if not issues else FAILURE_STATUS
    resolved_date_code = _extract_date_code(candidate_resolved, date_code)

    metrics: dict[str, Any] = {
        "status": status,
        "date_code": resolved_date_code,
        "reference": {
            "path": str(reference_resolved),
            "rows": int(len(reference_df)),
            "columns": list(reference_df.columns),
            "duplicates": reference_duplicates,
        },
        "candidate": {
            "path": str(candidate_resolved),
            "rows": int(len(candidate_df)),
            "columns": list(candidate_df.columns),
            "duplicates": candidate_duplicates,
        },
        "differences": {
            "cells_total": total_cells,
            "cells_different": cells_different,
            "rows_with_differences": rows_with_differences,
            "by_column": by_column,
        },
        "missing_rows": {
            "reference_only": missing_in_candidate,
            "candidate_only": missing_in_reference,
        },
        "issues": issues,
        "key_columns": list(key_columns),
    }

    report_stem = f"{REPORT_PREFIX}{resolved_date_code}"
    json_path = output_directory / f"{report_stem}{REPORT_JSON_SUFFIX}"
    markdown_path = output_directory / f"{report_stem}{REPORT_MD_SUFFIX}"

    with json_path.open("w", encoding=UTF8_ENCODING) as handle:
        json.dump(metrics, handle, indent=2, sort_keys=True)
        handle.write("\n")

    diff_path: Path | None = None
    if issues:
        diff_df = _differences_to_frame(
            mask,
            reference_aligned,
            candidate_aligned,
            key_columns=key_columns,
            limit=max_diff_rows,
        )
        if not diff_df.empty:
            diff_name = f"{DIFF_PREFIX}{resolved_date_code}{DIFF_SUFFIX}"
            diff_path = output_directory / diff_name
            write_csv_deterministic(
                diff_df,
                diff_path,
                encoding=UTF8_ENCODING,
                sep=delimiter,
                key_cols=list(key_columns),
            )
    markdown_path.write_text(
        _render_markdown(metrics, diff_path),
        encoding=UTF8_ENCODING,
    )

    return QaResult(
        passed=status == SUCCESS_STATUS,
        date_code=resolved_date_code,
        metrics=metrics,
        report_json=json_path,
        report_markdown=markdown_path,
        diff_csv=diff_path,
    )


# ===== CLI ==================================================================
def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate Python document post-processing against the reference output",
    )
    parser.add_argument(
        "--base-path",
        default=".",
        help="Root directory for resolving relative paths",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Relative or absolute path to the reference (Power Query) CSV",
    )
    parser.add_argument(
        "--actual",
        required=True,
        help="Relative or absolute path to the Python-generated CSV",
    )
    parser.add_argument(
        "--out",
        dest="output_dir",
        help="Optional directory for QA reports (defaults to the candidate directory)",
    )
    parser.add_argument(
        "--date-code",
        help="Override the date code extracted from the candidate filename",
    )
    parser.add_argument(
        "--delimiter",
        default=CSV_DELIMITER,
        help="Field delimiter used in the CSV files",
    )
    parser.add_argument(
        "--ref-encoding",
        default=CP1252_ENCODING,
        help="Encoding of the reference CSV",
    )
    parser.add_argument(
        "--actual-encoding",
        default=UTF8_ENCODING,
        help="Encoding of the Python-generated CSV",
    )
    parser.add_argument(
        "--max-diff-rows",
        type=int,
        default=DEFAULT_MAX_DIFF_ROWS,
        help="Number of discrepancies to include in the diff extract",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    result = run_document_postprocessing_check(
        base_path=args.base_path,
        reference_path=args.ref,
        candidate_path=args.actual,
        output_dir=args.output_dir,
        date_code=args.date_code,
        delimiter=args.delimiter,
        reference_encoding=args.ref_encoding,
        candidate_encoding=args.actual_encoding,
        max_diff_rows=args.max_diff_rows,
    )

    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
