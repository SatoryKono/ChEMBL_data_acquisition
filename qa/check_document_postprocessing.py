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
MAX_DIFF_KEY_EXPORT = 100
KEY_COLUMNS: tuple[str, str, str] = ("PMID", "document_chembl_id", "completed")
REPORT_PREFIX = "qa_document_postprocessing_report_"
DIFF_PREFIX = "qa_document_postprocessing_diff_"
REPORT_JSON_SUFFIX = ".json"
REPORT_MD_SUFFIX = ".md"
DIFF_SUFFIX = ".csv"
DATE_PATTERN = re.compile(r"(20\d{6})")
SUCCESS_STATUS = "passed"
FAILURE_STATUS = "failed"
BOOLEAN_TRUE_VALUES = frozenset({"true", "1", "yes"})
BOOLEAN_FALSE_VALUES = frozenset({"false", "0", "no", ""})
COMPLETED_FORMAT_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}$")
INVARIANT_SAMPLE_LIMIT = 5
INVARIANT_LABELS: Mapping[str, str] = {
    "invalid_formula": "invalid flag formula",
    "completed_format": "completed date format",
    "completed_order": "completed order monotonicity",
}


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
    key_count = 0
    for row_idx, row_mask in enumerate(mask.to_numpy()):
        if not row_mask.any():
            continue
        key_count += 1
        if key_count > limit:
            break
        index_value = mask.index[row_idx]
        if isinstance(index_value, tuple):
            key_map = dict(zip(key_columns, index_value, strict=False))
        else:
            key_map = {key_columns[0]: index_value}
            for extra in key_columns[1:]:
                key_map.setdefault(extra, "")
        for column_idx, has_difference in enumerate(row_mask):
            if not has_difference:
                continue
            record = {
                **key_map,
                "column": columns[column_idx],
                "expected": to_text(expected.iloc[row_idx, column_idx]),
                "actual": to_text(actual.iloc[row_idx, column_idx]),
            }
            records.append(record)
    return pd.DataFrame.from_records(records, columns=[*key_columns, "column", "expected", "actual"])


def _differences_by_column(mask: pd.DataFrame) -> dict[str, int]:
    """Return a dictionary counting mismatches per column."""

    return {column: int(mask[column].sum()) for column in mask.columns}


def _normalise_boolish(value: Any) -> bool | None:
    """Return a tri-state boolean for ``value`` supporting textual truthy forms."""

    text = to_text(value).strip().lower()
    if text in BOOLEAN_TRUE_VALUES:
        return True
    if text in BOOLEAN_FALSE_VALUES:
        return False
    return None


def _value_counts(series: pd.Series) -> Mapping[str, int]:
    """Return deterministic value counts for ``series`` as a mapping."""

    counts = series.fillna("").map(to_text).value_counts(sort=False)
    sorted_items = sorted(counts.items(), key=lambda item: item[0])
    return {str(key): int(value) for key, value in sorted_items}


def _summarise_invalid_flags(frame: pd.DataFrame) -> Mapping[str, Mapping[str, int]]:
    """Return per-column distributions for invalid flags within ``frame``."""

    columns = [column for column in frame.columns if column == "invalid" or column.startswith("invalid.")]
    summary: dict[str, Mapping[str, int]] = {}
    for column in sorted(columns):
        summary[column] = _value_counts(frame[column])
    return summary


def _boolean_share(series: pd.Series | None) -> Mapping[str, Any]:
    """Return counts and share of truthy values within ``series``."""

    if series is None:
        return {"available": False, "true": 0, "false": 0, "other": 0, "share_true": 0.0}

    values = series.fillna("").map(to_text)
    true_count = int(sum(1 for item in values if _normalise_boolish(item) is True))
    false_count = int(sum(1 for item in values if _normalise_boolish(item) is False))
    other_count = int(len(values) - true_count - false_count)
    total_rows = int(len(values))
    share_true = true_count / total_rows if total_rows else 0.0
    return {
        "available": True,
        "true": true_count,
        "false": false_count,
        "other": other_count,
        "share_true": share_true,
    }


def _collect_key_sample(index: pd.Index, key_columns: Sequence[str], limit: int = INVARIANT_SAMPLE_LIMIT) -> list[Mapping[str, Any]]:
    """Return a deterministic sample of key combinations from ``index``."""

    sample: list[Mapping[str, Any]] = []
    for value in index[:limit]:
        if isinstance(value, tuple):
            sample.append({column: element for column, element in zip(key_columns, value, strict=False)})
        else:
            row = {key_columns[0]: value}
            for extra in key_columns[1:]:
                row.setdefault(extra, "")
            sample.append(row)
    return sample


def _check_invalid_formula(frame: pd.DataFrame, key_columns: Sequence[str]) -> Mapping[str, Any]:
    """Return invariant status ensuring ``invalid`` equals the OR of component flags."""

    if "invalid" not in frame.columns:
        return {"available": False, "passed": True, "violations": 0, "failing_keys": []}

    detail_columns = [column for column in frame.columns if column.startswith("invalid.")]
    invalid_true = frame["invalid"].map(_normalise_boolish).eq(True)
    if detail_columns:
        detail_true = frame[detail_columns].map(_normalise_boolish).eq(True)
        expected_true = detail_true.any(axis=1)
    else:
        expected_true = pd.Series(False, index=frame.index)
    violations_mask = invalid_true != expected_true
    violations = int(violations_mask.sum())
    failing_keys = _collect_key_sample(frame.index[violations_mask], key_columns)
    return {
        "available": True,
        "passed": violations == 0,
        "violations": violations,
        "failing_keys": failing_keys,
    }


def _check_completed_format(frame: pd.DataFrame, key_columns: Sequence[str]) -> Mapping[str, Any]:
    """Return invariant status confirming ``completed`` follows the canonical pattern."""

    if "completed" not in frame.columns:
        return {"available": False, "passed": True, "violations": 0, "failing_keys": []}

    completed_values = frame["completed"].map(to_text)
    valid_mask = completed_values.str.fullmatch(COMPLETED_FORMAT_PATTERN) | completed_values.eq("")
    violations_mask = ~valid_mask
    violations = int(violations_mask.sum())
    failing_keys = _collect_key_sample(frame.index[violations_mask], key_columns)
    return {
        "available": True,
        "passed": violations == 0,
        "violations": violations,
        "failing_keys": failing_keys,
    }


def _parse_completed_tuple(value: str) -> tuple[int, int, int] | None:
    """Return a sortable tuple extracted from ``value`` when formatted correctly."""

    if not COMPLETED_FORMAT_PATTERN.fullmatch(value):
        return None
    parts = value.split("-")
    return tuple(int(part) for part in parts)


def _check_completed_order(frame: pd.DataFrame, key_columns: Sequence[str]) -> Mapping[str, Any]:
    """Return invariant status ensuring completed dates are non-decreasing by sort order."""

    if "completed" not in frame.columns or "sortorder.document" not in frame.columns:
        return {"available": False, "passed": True, "violations": 0, "failing_keys": []}

    ordered = frame[["completed", "sortorder.document"]].copy()
    ordered["completed"] = ordered["completed"].map(to_text)
    ordered["sortorder.document"] = ordered["sortorder.document"].map(to_text)
    ordered = ordered.sort_values("sortorder.document", kind="stable")

    previous: tuple[int, int, int] | None = None
    violating_labels: list[Any] = []
    for label, value in ordered["completed"].items():
        parsed = _parse_completed_tuple(value)
        if parsed is None:
            continue
        if previous is None:
            previous = parsed
            continue
        if parsed < previous:
            violating_labels.append(label)
        else:
            previous = parsed
    violations = len(violating_labels)
    failing_keys = _collect_key_sample(pd.Index(violating_labels), key_columns)
    return {
        "available": True,
        "passed": violations == 0,
        "violations": violations,
        "failing_keys": failing_keys,
    }


def _summarise_dataset(
    *,
    frame: pd.DataFrame,
    canonical: pd.DataFrame,
    key_columns: Sequence[str],
    path: str,
    duplicate_count: int,
) -> Mapping[str, Any]:
    """Return derived metrics for ``frame`` and ``canonical`` representations."""

    review_summary = _boolean_share(canonical.get("review"))
    experimental_summary = _boolean_share(canonical.get("experimental"))
    invariants = {
        "invalid_formula": _check_invalid_formula(canonical, key_columns),
        "completed_format": _check_completed_format(canonical, key_columns),
        "completed_order": _check_completed_order(canonical, key_columns),
    }
    return {
        "path": path,
        "rows": int(len(frame)),
        "columns": list(frame.columns),
        "column_count": int(len(frame.columns)),
        "duplicates": duplicate_count,
        "invalid_flags": _summarise_invalid_flags(canonical),
        "review": review_summary,
        "experimental": experimental_summary,
        "invariants": invariants,
    }


def _render_markdown(metrics: Mapping[str, Any], diff_path: Path | None) -> str:
    """Return a Markdown summary of the QA run."""

    def _format_bool(flag: bool) -> str:
        return "✅ yes" if flag else "❌ no"

    def _format_invalid(summary: Mapping[str, Mapping[str, int]]) -> str:
        if not summary:
            return "n/a"
        parts = []
        for column, values in summary.items():
            value_text = ", ".join(
                f"{name or '""'}={count}" for name, count in values.items()
            )
            parts.append(f"{column} [{value_text}]")
        return "; ".join(parts)

    def _format_share_block(summary: Mapping[str, Any], total_rows: int) -> str:
        if not summary.get("available"):
            return "n/a"
        true_count = summary["true"]
        percentage = summary["share_true"] * 100
        return f"{percentage:.2f}% ({true_count}/{total_rows})"

    ref_rows = metrics["reference"]["rows"]
    cand_rows = metrics["candidate"]["rows"]
    lines = ["# Document post-processing QA report", ""]
    lines.append(f"- Status: **{metrics['status']}**")
    lines.append(f"- Reference rows: {ref_rows}")
    lines.append(f"- Candidate rows: {cand_rows}")
    lines.append(
        "- Column sets identical: "
        + _format_bool(metrics["structure"]["columns_equal"])
    )
    lines.append(
        "- Column order identical: "
        + _format_bool(metrics["structure"]["column_order_equal"])
    )
    lines.append(
        "- Cells different: "
        f"{metrics['differences']['cells_different']} / {metrics['differences']['cells_total']} "
        f"(rows impacted: {metrics['differences']['rows_with_differences']} of {metrics['differences']['rows_compared']})"
    )
    lines.append(
        "- Missing keys: "
        f"reference-only {metrics['missing_rows']['reference_only']['count']}, "
        f"candidate-only {metrics['missing_rows']['candidate_only']['count']}"
    )
    lines.append(
        "- Invalid flags (reference): "
        + _format_invalid(metrics["reference"]["invalid_flags"])
    )
    lines.append(
        "- Invalid flags (candidate): "
        + _format_invalid(metrics["candidate"]["invalid_flags"])
    )
    lines.append(
        "- Review share: reference "
        + _format_share_block(metrics["reference"]["review"], ref_rows)
        + "; candidate "
        + _format_share_block(metrics["candidate"]["review"], cand_rows)
    )
    lines.append(
        "- Experimental share: reference "
        + _format_share_block(metrics["reference"]["experimental"], ref_rows)
        + "; candidate "
        + _format_share_block(metrics["candidate"]["experimental"], cand_rows)
    )
    invariants = metrics["candidate"]["invariants"]
    ref_invariants = metrics["reference"]["invariants"]
    lines.append(
        "- Invariants: "
        f"invalid formula ref {_format_bool(ref_invariants['invalid_formula']['passed'])} / "
        f"cand {_format_bool(invariants['invalid_formula']['passed'])}; "
        f"completed format ref {_format_bool(ref_invariants['completed_format']['passed'])} / "
        f"cand {_format_bool(invariants['completed_format']['passed'])}; "
        f"order ref {_format_bool(ref_invariants['completed_order']['passed'])} / "
        f"cand {_format_bool(invariants['completed_order']['passed'])}"
    )
    lines.append(f"- Issues detected: {len(metrics['issues'])}")
    if diff_path is not None:
        lines.append(
            f"- Diff excerpt: [{diff_path.name}]({diff_path.name})"
        )
    else:
        lines.append("- Diff excerpt: not generated")

    if metrics["issues"]:
        lines.append("")
        lines.append("## Issues")
        for issue in metrics["issues"]:
            lines.append(f"- {issue}")

    if metrics["differences"]["by_column"]:
        lines.append("")
        lines.append("## Differences by column")
        lines.append("| Column | Differences |")
        lines.append("| --- | ---: |")
        for column, count in metrics["differences"]["by_column"].items():
            lines.append(f"| {column} | {count} |")

    lines.append("")
    return "\n".join(lines)


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

    structure_metrics = {
        "columns_equal": set(reference_df.columns) == set(candidate_df.columns),
        "column_order_equal": list(reference_df.columns) == list(candidate_df.columns),
    }

    reference_summary = _summarise_dataset(
        frame=reference_df,
        canonical=reference_canonical,
        key_columns=key_columns,
        path=str(reference_resolved),
        duplicate_count=reference_duplicates,
    )
    candidate_summary = _summarise_dataset(
        frame=candidate_df,
        canonical=candidate_canonical,
        key_columns=key_columns,
        path=str(candidate_resolved),
        duplicate_count=candidate_duplicates,
    )

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
    if structure_metrics["columns_equal"] and not structure_metrics["column_order_equal"]:
        issues.append("Column order differs between reference and candidate outputs")

    for dataset_name, dataset_summary in (
        ("reference", reference_summary),
        ("candidate", candidate_summary),
    ):
        for invariant_key, invariant in dataset_summary["invariants"].items():
            if invariant.get("available") and not invariant.get("passed"):
                label = INVARIANT_LABELS.get(invariant_key, invariant_key.replace("_", " "))
                issues.append(
                    f"{dataset_name.title()} {label} violated for {invariant['violations']} rows"
                )

    status = SUCCESS_STATUS if not issues else FAILURE_STATUS
    resolved_date_code = _extract_date_code(candidate_resolved, date_code)

    metrics: dict[str, Any] = {
        "status": status,
        "date_code": resolved_date_code,
        "reference": reference_summary,
        "candidate": candidate_summary,
        "structure": structure_metrics,
        "differences": {
            "cells_total": total_cells,
            "cells_different": cells_different,
            "rows_compared": int(len(mask.index)),
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
        diff_limit = min(max_diff_rows, MAX_DIFF_KEY_EXPORT)
        diff_df = _differences_to_frame(
            mask,
            reference_aligned,
            candidate_aligned,
            key_columns=key_columns,
            limit=diff_limit,
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
        help="Root directory containing the Power Query reference and Python outputs",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Relative or absolute path to the Python-generated document CSV (output.document_*.csv)",
    )
    parser.add_argument(
        "--reports-dir",
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

    base_dir = Path(args.base_path).resolve()
    reference_path = _resolve_relative(base_dir, Path("input") / "full" / "document.csv")
    candidate_path = _resolve_relative(base_dir, args.out)

    result = run_document_postprocessing_check(
        base_path=base_dir,
        reference_path=reference_path,
        candidate_path=candidate_path,
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
