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
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from library import document_postprocessing as dp
from library.config import IoCfg
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
    """Compare the Power Query reproduction with the Python post-processing result."""

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

    ref_source = _load_csv(
        reference_resolved,
        encoding=reference_encoding,
        delimiter=delimiter,
        keep_default_na=True,
    )
    out_source = _load_csv(
        candidate_resolved,
        encoding=candidate_encoding,
        delimiter=delimiter,
        keep_default_na=True,
    )

    expected_df = _power_query_expected(ref_source, out_source)

    cfg = IoCfg(csv_sep=delimiter, csv_encoding=candidate_encoding)
    processed_path = dp.postprocess_file(
        candidate_resolved,
        candidate_resolved.parent,
        cfg=cfg,
        ref_document_path=reference_resolved,
    )

    actual_df = _load_csv(
        processed_path,
        encoding=dp.UTF8_ENCODING,
        delimiter=delimiter,
        keep_default_na=False,
    )

    expected_df = _ensure_columns(expected_df, key_columns)
    actual_df = _ensure_columns(actual_df, key_columns)

    reference_canonical = _canonicalise(expected_df)
    candidate_canonical = _canonicalise(actual_df)

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

    missing_in_actual = sorted(set(all_columns) - set(actual_df.columns))
    missing_in_expected = sorted(set(all_columns) - set(expected_df.columns))

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
    if missing_in_expected:
        issues.append(
            "Reference missing columns present in candidate: "
            + ", ".join(missing_in_expected)
        )
    if missing_in_actual:
        issues.append(
            "Candidate missing columns present in reference: "
            + ", ".join(missing_in_actual)
        )

    status = SUCCESS_STATUS if not issues else FAILURE_STATUS
    resolved_date_code = _extract_date_code(processed_path, date_code)

    metrics: dict[str, Any] = {
        "status": status,
        "date_code": resolved_date_code,
        "reference": {
            "path": str(reference_resolved),
            "rows": int(len(expected_df)),
            "columns": list(expected_df.columns),
            "duplicates": reference_duplicates,
        },
        "candidate": {
            "path": str(processed_path),
            "rows": int(len(actual_df)),
            "columns": list(actual_df.columns),
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
        help="Relative or absolute path to the unprocessed out_document CSV",
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
# ===== Helpers ===============================================================
def _load_csv(
    path: Path,
    *,
    encoding: str,
    delimiter: str,
    keep_default_na: bool,
) -> pd.DataFrame:
    """Return a DataFrame loaded from ``path`` with consistent options."""

    return pd.read_csv(
        path,
        dtype="string",
        encoding=encoding,
        sep=delimiter,
        keep_default_na=keep_default_na,
    )


def _prepare_reference_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` shaped like the Power Query ``ref_document`` table."""

    reference = frame.copy()

    if "classification" in reference.columns:
        classification = (
            pd.to_numeric(reference["classification"], errors="coerce")
            .fillna(0)
            .astype(int)
        )
        reference["doctype_review"] = classification.astype(bool)
    elif "doctype_review" not in reference.columns:
        reference["doctype_review"] = False

    drop_columns = [
        "abstract",
        "authors",
        "DOI",
        "first_page",
        "issue",
        "journal",
        "last_page",
        "month",
        "postcodes",
        "title",
        "volume",
        "year",
        "classification",
    ]
    reference = reference.drop(columns=[c for c in drop_columns if c in reference.columns])

    reference["document_chembl_id"] = reference["document_chembl_id"].map(to_text)

    return reference


def _series_from_column(
    frame: pd.DataFrame,
    column: str,
    *,
    transform: Callable[[Any], Any] | None = to_text,
    default: str = "",
) -> pd.Series:
    """Return ``column`` mapped through ``transform`` with ``default`` fallback."""

    if column in frame.columns:
        series = frame[column]
    else:
        series = pd.Series([default] * len(frame), index=frame.index)
    if transform is None:
        return series
    return series.map(transform)


def _power_query_expected(
    ref_document: pd.DataFrame,
    out_document: pd.DataFrame,
) -> pd.DataFrame:
    """Return the expected Power Query output for ``out_document``."""

    if out_document.empty:
        return pd.DataFrame(columns=list(dp.FINAL_COLUMN_ORDER))

    frame = out_document.copy()

    for column in frame.columns:
        frame[column] = frame[column].map(to_text)

    for column in dp.PUBMED_NUMERIC_TEXT_COLUMNS:
        frame[column] = _series_from_column(frame, column).replace("", "0")

    for column in dp.CHEMBL_NUMERIC_TEXT_COLUMNS + dp.PUBMED_NUMERIC_TEXT_COLUMNS:
        frame[column] = _series_from_column(frame, column)

    frame["ChEMBL.journal"] = _series_from_column(
        frame, "ChEMBL.journal", transform=dp.normalize_journal
    )
    frame["PubMed.JournalTitle"] = _series_from_column(
        frame, "PubMed.JournalTitle", transform=dp.normalize_journal
    )

    ref_doi = _series_from_column(frame, "ChEMBL.doi")
    doi_agree: list[pd.Series] = []
    doi_conflict: list[pd.Series] = []
    for column in dp.EXTERNAL_DOI_COLUMNS:
        external = _series_from_column(frame, column)
        valid = ~external.isin(("", "0"))
        doi_agree.append(valid & (external == ref_doi))
        doi_conflict.append(valid & (external != ref_doi))
    agree_any = dp._series_any(doi_agree) if doi_agree else pd.Series(False, index=frame.index)
    conflict_any = (
        dp._series_any(doi_conflict)
        if doi_conflict
        else pd.Series(False, index=frame.index)
    )
    frame["invalid.doi"] = (~agree_any) & conflict_any

    ref_pmid_numbers = _series_from_column(frame, "ChEMBL.pubmed_id").map(dp.try_parse_int)
    ref_valid = ref_pmid_numbers.notna()
    pmid_agree: list[pd.Series] = []
    pmid_conflict: list[pd.Series] = []
    for column in dp.EXTERNAL_PMID_COLUMNS:
        external_text = _series_from_column(frame, column)
        external_numbers = external_text.map(dp.try_parse_int)
        external_valid = external_numbers.notna()
        match = ref_valid & external_valid & (external_numbers == ref_pmid_numbers)
        pmid_agree.append(match)

        non_empty = external_text != ""
        mismatch = (
            non_empty & ref_valid & external_valid & (external_numbers != ref_pmid_numbers)
        )
        pmid_conflict.append(mismatch)
    agree_pmid_any = (
        dp._series_any(pmid_agree)
        if pmid_agree
        else pd.Series(False, index=frame.index)
    )
    conflict_pmid_any = (
        dp._series_any(pmid_conflict)
        if pmid_conflict
        else pd.Series(False, index=frame.index)
    )
    frame["invalid.PMID"] = (~agree_pmid_any) & conflict_pmid_any

    review_columns = (
        "PubMed.PublicationType",
        "scholar.PublicationTypes",
        "OpenAlex.PublicationTypes",
        "OpenAlex.TypeCrossref",
        "crossref.Type",
    )
    review_series = [
        _series_from_column(frame, column)
        .astype("string")
        .fillna("")
        .str.lower()
        .str.contains("review", na=False)
        for column in review_columns
    ]
    frame["review"] = dp._series_any(review_series) if review_series else pd.Series(False, index=frame.index)

    journal_match = (
        frame["PubMed.JournalTitle"] == frame["ChEMBL.journal"]
    ) & (frame["ChEMBL.journal"] != "")
    volume_match = (
        frame["PubMed.Volume"] == frame["ChEMBL.volume"]
    ) & (frame["ChEMBL.volume"] != "")
    issue_match = (
        frame["PubMed.Issue"] == frame["ChEMBL.issue"]
    ) & (frame["ChEMBL.issue"] != "")
    start_match = (
        frame["PubMed.StartPage"] == frame["ChEMBL.first_page"]
    ) & (frame["ChEMBL.first_page"] != "")
    end_match = (
        frame["PubMed.EndPage"] == frame["ChEMBL.last_page"]
    ) & (frame["ChEMBL.last_page"] != "")

    journal_value = np.where(journal_match, frame["ChEMBL.journal"], "unknown")
    volume_value = np.where(volume_match, frame["ChEMBL.volume"], "unknown")
    issue_value = np.where(issue_match, frame["ChEMBL.issue"], "unknown")
    start_value = np.where(start_match, frame["ChEMBL.first_page"], "unknown")
    end_value = np.where(end_match, frame["ChEMBL.last_page"], "unknown")

    frame["reference"] = (
        pd.Series(journal_value, index=frame.index)
        + ", "
        + pd.Series(volume_value, index=frame.index)
        + "("
        + pd.Series(issue_value, index=frame.index)
        + "), p."
        + pd.Series(start_value, index=frame.index)
        + "-"
        + pd.Series(end_value, index=frame.index)
    )

    frame["invalid.reference"] = frame["reference"].str.contains("unknown", na=False)

    year_completed = _series_from_column(frame, "PubMed.YearCompleted")
    year_revised = _series_from_column(frame, "PubMed.YearRevised")
    chembl_year = _series_from_column(frame, "ChEMBL.year")
    frame["year"] = np.where(
        ~year_completed.map(dp.null_or_empty),
        year_completed.map(dp.pad4),
        np.where(
            ~year_revised.map(dp.null_or_empty),
            year_revised.map(dp.pad4),
            chembl_year.map(dp.pad4),
        ),
    )

    month_completed = _series_from_column(frame, "PubMed.MonthCompleted")
    month_revised = _series_from_column(frame, "PubMed.MonthRevised")
    frame["month"] = np.where(
        ~month_completed.map(dp.null_or_empty),
        month_completed.map(dp.pad2),
        np.where(
            ~month_revised.map(dp.null_or_empty),
            month_revised.map(dp.pad2),
            "00",
        ),
    )

    day_completed = _series_from_column(frame, "PubMed.DayCompleted")
    day_revised = _series_from_column(frame, "PubMed.DayRevised")
    frame["day"] = np.where(
        ~day_completed.map(dp.null_or_empty),
        day_completed.map(dp.pad2),
        np.where(
            ~day_revised.map(dp.null_or_empty),
            day_revised.map(dp.pad2),
            "00",
        ),
    )

    frame["completed"] = (
        pd.Series(frame["year"], index=frame.index)
        + "-"
        + pd.Series(frame["month"], index=frame.index)
        + "-"
        + pd.Series(frame["day"], index=frame.index)
    )

    frame["ChEMBL.pubmed_id"] = _series_from_column(frame, "ChEMBL.pubmed_id").map(
        dp.pad_pmid8
    )
    frame["sortorder.document"] = (
        _series_from_column(frame, "PubMed.ISSN")
        + ":"
        + frame["completed"]
        + ":"
        + frame["ChEMBL.pubmed_id"].map(to_text)
    )

    frame = frame.rename(columns={"PubMed.PMID": "PMID", "ChEMBL.doi": "doi"})
    frame["PMID"] = frame["PMID"].map(to_text)
    frame["doi"] = frame["doi"].map(to_text)

    frame["invalid"] = (
        frame["invalid.doi"] | frame["invalid.PMID"] | frame["invalid.reference"]
    )

    drop_candidates = [
        column for column in dp.STAGE_REMOVED_COLUMNS if column in frame.columns
    ]
    frame = frame.drop(columns=drop_candidates)

    rename_map = {
        "ChEMBL.document_chembl_id": "document_chembl_id",
        "PubMed.DOI": "doi",
        "PubMed.ArticleTitle": "title",
        "PubMed.Abstract": "abstract",
        "PubMed.MeSH_Descriptors": "mesh.descriptors",
        "PubMed.MeSH_Qualifiers": "mesh.qualifiers",
        "PubMed.ChemicalList": "chemical_list",
        "OpenAlex.MeshDescriptors": "OpenAlex.mesh.descriptors",
    }
    frame = frame.rename(columns=rename_map)
    frame["doi"] = frame["doi"].map(to_text)

    secondary_rename = {
        "chemical_list": "PubMed.chemical.list",
        "mesh.qualifiers": "PubMed.mesh.qualifiers",
        "mesh.descriptors": "PubMed.mesh.descriptors",
    }
    frame = frame.rename(columns=secondary_rename)

    for column in dp.LOWERCASE_LIST_COLUMNS:
        if column in frame.columns:
            frame[column] = frame[column].map(dp.safe_lower)

    harmonised_reference = _prepare_reference_frame(ref_document)

    frame = frame.merge(harmonised_reference, on="document_chembl_id", how="left")

    if "document_contains_external_links" in frame.columns:
        def _normalize_bool(value: Any) -> Any:
            if pd.isna(value):
                return pd.NA
            text = to_text(value).strip().lower()
            if text in {"", "nan"}:
                return pd.NA
            if text in {"true", "1", "yes"}:
                return True
            if text in {"false", "0", "no"}:
                return False
            return bool(value)

        frame["document_contains_external_links"] = frame[
            "document_contains_external_links"
        ].map(_normalize_bool)

    review_values: list[Any] = []
    for current, doctype in zip(frame["review"], frame["doctype_review"]):
        current_value = None if pd.isna(current) else bool(current)
        doctype_value = None if pd.isna(doctype) else bool(doctype)
        if current_value is True or doctype_value is True:
            review_values.append(True)
        elif current_value is False and doctype_value is False:
            review_values.append(False)
        else:
            review_values.append(pd.NA)
    frame["review"] = pd.Series(review_values, index=frame.index, dtype="boolean")

    experimental_values: list[Any] = []
    for value in frame["review"]:
        if value is pd.NA or pd.isna(value):
            experimental_values.append(pd.NA)
        else:
            experimental_values.append(not bool(value))
    frame["experimental"] = pd.Series(
        experimental_values, index=frame.index, dtype="boolean"
    )

    frame = frame.drop(columns=["is_experimental_doc"], errors="ignore")

    boolean_columns = [
        "review",
        "experimental",
        "invalid",
        "invalid.doi",
        "invalid.PMID",
        "invalid.reference",
        "document_contains_external_links",
    ]
    for column in boolean_columns:
        if column in frame.columns:
            frame[column] = frame[column].map(
                lambda value: ""
                if pd.isna(value)
                else str(bool(value)).lower()
            )

    missing = [column for column in dp.FINAL_COLUMN_ORDER if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing expected columns after harmonisation: {missing}")

    frame = frame.loc[:, dp.FINAL_COLUMN_ORDER]
    frame = frame.sort_values("completed", ascending=True, kind="mergesort")
    frame.reset_index(drop=True, inplace=True)

    return frame
