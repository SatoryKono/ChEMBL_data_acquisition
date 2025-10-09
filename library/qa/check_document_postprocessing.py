"""Crosswalk-based QA for document post-processing outputs.

Changelog
~~~~~~~~
- 2025-02-14: Replaced legacy CSV diff with crosswalk-driven comparison of the
  Python projection against the canonical Power Query logic.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import pandas as pd
import yaml

from library.common.csv_utils import write_csv_deterministic
from library.common.text_utils import to_text


# ===== Parameters ============================================================
CP1252_ENCODING = "cp1252"
UTF8_ENCODING = "utf-8"
CSV_DELIMITER = ","
DEFAULT_BASE_PATH = Path("data")
DEFAULT_REPORT_DIR = Path("output") / "document"
CROSSWALK_PATH = Path(__file__).resolve().parent / "document_schema_crosswalk.yaml"
DATE_PATTERN = re.compile(r"(20\d{6})")
DEFAULT_DIFF_LIMIT = 100
MAX_DIFF_KEY_EXPORT = DEFAULT_DIFF_LIMIT
BOOLEAN_TRUE_VALUES = {"true", "1", "yes", "y", "t"}
BOOLEAN_FALSE_VALUES = {"false", "0", "no", "n", "f", ""}
MESH_SPLIT_PATTERN = re.compile(r"[;,|]")
PROVIDER_COLUMNS_M: Mapping[str, tuple[str, ...]] = {
    "chembl": (
        "ChEMBL.title",
        "ChEMBL.abstract",
        "ChEMBL.doi",
        "ChEMBL.year",
        "ChEMBL.journal",
        "ChEMBL.journal_abbrev",
        "ChEMBL.volume",
        "ChEMBL.issue",
        "ChEMBL.first_page",
        "ChEMBL.last_page",
        "ChEMBL.pubmed_id",
        "ChEMBL.authors",
        "ChEMBL.source",
        "title",
        "abstract",
        "doi",
    ),
    "pubmed": (
        "PMID",
        "PubMed.PMID",
        "PubMed.ArticleTitle",
        "PubMed.Abstract",
        "PubMed.mesh.descriptors",
        "PubMed.mesh.qualifiers",
        "PubMed.chemical.list",
        "PubMed.Volume",
        "PubMed.Issue",
        "PubMed.StartPage",
        "PubMed.EndPage",
        "PubMed.ISSN",
        "PubMed.PublicationType",
        "PubMed.YearCompleted",
        "PubMed.MonthCompleted",
        "PubMed.DayCompleted",
        "PubMed.YearRevised",
        "PubMed.MonthRevised",
        "PubMed.DayRevised",
        "PubMed.Error",
    ),
    "semantic_scholar": (
        "scholar.PMID",
        "scholar.DOI",
        "scholar.PublicationTypes",
        "scholar.Venue",
        "scholar.SemanticScholarId",
        "scholar.ExternalIds",
        "scholar.Error",
    ),
    "openalex": (
        "OpenAlex.PMID",
        "OpenAlex.DOI",
        "OpenAlex.TypeCrossref",
        "OpenAlex.Id",
        "OpenAlex.mesh.descriptors",
        "OpenAlex.mesh.qualifiers",
        "OpenAlex.Error",
    ),
    "crossref": (
        "crossref.DOI",
        "crossref.Type",
        "crossref.Title",
        "crossref.Error",
    ),
}
INVALID_ERROR_TOKEN_MAP = {
    "invalid.doi": {"chembl", "pubmed", "semantic_scholar", "openalex", "crossref", "unknown", "doi"},
    "invalid.PMID": {"pubmed", "semantic_scholar", "openalex", "unknown", "pmid"},
}
MATCH_THRESHOLDS = {
    "ids": 0.995,
    "preferred": 0.995,
    "publication_year": 0.995,
    "is_review": 1.0,
    "mesh_terms": 0.995,
    "has_flags": 0.995,
    "coverage": 0.995,
}
REPORT_PREFIX = "qa_document_postprocessing_report_"
REPORT_JSON_SUFFIX = ".json"
REPORT_MD_SUFFIX = ".md"
DIFF_PREFIX = "qa_document_postprocessing_diff_"
DIFF_SUFFIX = ".csv"


# ===== Crosswalk structures ==================================================
@dataclass(frozen=True)
class BuilderConfig:
    """Configuration for deriving a comparable series."""

    builder: str
    columns: tuple[str, ...] = field(default_factory=tuple)
    params: Mapping[str, Any] = field(default_factory=dict)

    @staticmethod
    def from_mapping(data: Mapping[str, Any]) -> "BuilderConfig":
        builder = data.get("builder")
        if not builder:
            raise ValueError("builder key is required in crosswalk field")
        columns = tuple(data.get("columns", ()))
        params = data.get("params", {})
        return BuilderConfig(builder=builder, columns=columns, params=params)


@dataclass(frozen=True)
class CrosswalkField:
    """Crosswalk definition for a single comparable field."""

    name: str
    description: str
    metric_group: str
    python: BuilderConfig
    m: BuilderConfig

    @staticmethod
    def from_mapping(data: Mapping[str, Any]) -> "CrosswalkField":
        name = data.get("name")
        if not name:
            raise ValueError("Field name missing in crosswalk")
        description = data.get("description", "")
        metric_group = data.get("metric_group", "misc")
        python_cfg = BuilderConfig.from_mapping(data.get("python", {}))
        m_cfg = BuilderConfig.from_mapping(data.get("m", {}))
        return CrosswalkField(
            name=name,
            description=description,
            metric_group=metric_group,
            python=python_cfg,
            m=m_cfg,
        )


@dataclass(frozen=True)
class Crosswalk:
    """Container for the full crosswalk specification."""

    version: str
    description: str
    fields: tuple[CrosswalkField, ...]

    @staticmethod
    def load(path: Path) -> "Crosswalk":
        with path.open("r", encoding=UTF8_ENCODING) as stream:
            payload = yaml.safe_load(stream)
        version = payload.get("version", "unknown")
        description = payload.get("description", "")
        fields = tuple(CrosswalkField.from_mapping(item) for item in payload.get("fields", []))
        if not fields:
            raise ValueError("Crosswalk must define at least one field")
        return Crosswalk(version=version, description=description, fields=fields)

    def metric_groups(self) -> Mapping[str, list[str]]:
        groups: dict[str, list[str]] = {}
        for field_cfg in self.fields:
            groups.setdefault(field_cfg.metric_group, []).append(field_cfg.name)
        return groups


# ===== Helpers ===============================================================
def _resolve_path(base: Path, candidate: str | os.PathLike[str]) -> Path:
    text = str(candidate)
    candidate_path = Path(text.replace("\\", os.sep))
    if candidate_path.is_absolute():
        return candidate_path
    return (base / candidate_path).resolve()


def _read_csv(path: Path, *, encoding: str) -> pd.DataFrame:
    return pd.read_csv(
        path,
        sep=CSV_DELIMITER,
        encoding=encoding,
        dtype=str,
        keep_default_na=False,
    )


def _string_series(series: pd.Series) -> pd.Series:
    if series is None:
        return pd.Series(dtype="string")
    return series.fillna("").astype("string").str.strip()


def _normalize_numeric_identifier(value: Any) -> str:
    text = to_text(value).strip()
    if not text:
        return ""
    digits = "".join(ch for ch in text if ch.isdigit())
    if not digits:
        return ""
    return str(int(digits))


def _normalize_numeric_series(series: pd.Series) -> pd.Series:
    return series.map(_normalize_numeric_identifier)


def _normalize_boolean_value(value: Any) -> bool | None:
    text = to_text(value).strip().lower()
    if text in BOOLEAN_TRUE_VALUES:
        return True
    if text in BOOLEAN_FALSE_VALUES:
        return False
    return None


def _boolean_series(series: pd.Series) -> pd.Series:
    normalized = series.map(_normalize_boolean_value)
    return normalized.astype("boolean")


def _mesh_tokens(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    tokens: list[str] = []
    for value in values:
        for token in MESH_SPLIT_PATTERN.split(value):
            cleaned = token.strip().lower()
            if not cleaned:
                continue
            if cleaned in seen:
                continue
            seen.add(cleaned)
            tokens.append(cleaned)
    tokens.sort()
    return tokens


def _mesh_series(columns: Sequence[pd.Series]) -> pd.Series:
    if not columns:
        return pd.Series(dtype="string")
    length = len(columns[0])
    result = []
    for row_idx in range(length):
        values = [to_text(column.iloc[row_idx]) for column in columns]
        tokens = _mesh_tokens(values)
        result.append("; ".join(tokens))
    return pd.Series(result, dtype="string")


def _year_from_date(value: str) -> str:
    text = to_text(value).strip()
    if not text:
        return ""
    match = re.match(r"(\d{4})", text)
    return match.group(1) if match else ""


def _coalesce_columns(frame: pd.DataFrame, columns: Sequence[str]) -> pd.Series:
    if not columns:
        return pd.Series([""] * len(frame), index=frame.index, dtype="string")
    result = pd.Series([""] * len(frame), index=frame.index, dtype="string")
    for column in columns:
        if column not in frame.columns:
            continue
        candidate = _string_series(frame[column])
        mask = result.eq("") & candidate.ne("")
        result.loc[mask] = candidate.loc[mask]
        if result.ne("").all():
            break
    return result


def _tokenise_sources(value: str) -> set[str]:
    tokens = MESH_SPLIT_PATTERN.split(to_text(value).lower())
    return {token.strip() for token in tokens if token.strip()}


def _extract_date_code(path: Path) -> str:
    match = DATE_PATTERN.search(path.name)
    if match:
        return match.group(1)
    return datetime.utcnow().strftime("%Y%m%d")


def _collect_key_sample(index: pd.Index, limit: int = 5) -> list[Mapping[str, str]]:
    sample: list[Mapping[str, str]] = []
    for key in index[:limit]:
        if isinstance(key, tuple):
            sample.append({"document_chembl_id": to_text(key[0]), "primary_pubmed_id": to_text(key[1])})
        else:
            sample.append({"document_chembl_id": to_text(key), "primary_pubmed_id": ""})
    return sample


# ===== Dataset contexts ======================================================
class DatasetContext:
    """Base class exposing helper utilities for builder functions."""

    def __init__(self, frame: pd.DataFrame):
        self.frame = frame

    @property
    def length(self) -> int:
        return len(self.frame)


class PythonDatasetContext(DatasetContext):
    """Context for the Python post-processed dataset."""

    def get_column(self, column: str) -> pd.Series:
        if column not in self.frame.columns:
            return pd.Series(["" for _ in range(self.length)], dtype="string")
        return _string_series(self.frame[column])

    def provider_flag(self, provider: str) -> pd.Series:
        column_name = f"has_{provider}"
        series = self.frame.get(column_name)
        if series is None:
            return pd.Series([False] * self.length, dtype="boolean")
        values = _boolean_series(_string_series(series))
        return values.fillna(False)

    def metadata_count(self) -> pd.Series:
        series = self.frame.get("metadata_source_count")
        if series is None:
            return pd.Series([0] * self.length, dtype="Int64")
        numeric = pd.to_numeric(_string_series(series).replace("", pd.NA), errors="coerce")
        return numeric.astype("Int64").fillna(0)


class MPowerQueryContext(DatasetContext):
    """Context for the legacy Power Query shaped dataset."""

    def __init__(self, frame: pd.DataFrame):
        super().__init__(frame)
        self._provider_cache: dict[str, pd.Series] = {}

    def get_column(self, column: str) -> pd.Series:
        if column not in self.frame.columns:
            return pd.Series(["" for _ in range(self.length)], dtype="string")
        return _string_series(self.frame[column])

    def provider_flag(self, provider: str) -> pd.Series:
        if provider in self._provider_cache:
            return self._provider_cache[provider]
        columns = PROVIDER_COLUMNS_M.get(provider, ())
        flag = pd.Series([False] * self.length, dtype="boolean")
        for column in columns:
            if column not in self.frame.columns:
                continue
            candidate = _string_series(self.frame[column]).ne("")
            flag = flag | candidate
        self._provider_cache[provider] = flag.fillna(False)
        return self._provider_cache[provider]

    def metadata_count(self) -> pd.Series:
        flags = [self.provider_flag(provider) for provider in PROVIDER_COLUMNS_M]
        if not flags:
            return pd.Series([0] * self.length, dtype="Int64")
        total = pd.DataFrame(flags).T.sum(axis=1)
        return total.astype("Int64")


# ===== Builder registries ====================================================
def _build_python_series(context: PythonDatasetContext, config: BuilderConfig) -> pd.Series:
    if config.builder == "column":
        return context.get_column(config.columns[0]) if config.columns else context.get_column("")
    if config.builder == "numeric_identifier":
        if not config.columns:
            return pd.Series(["" for _ in range(context.length)], dtype="string")
        series = context.get_column(config.columns[0])
        return _normalize_numeric_series(series)
    if config.builder == "boolean":
        if not config.columns:
            return pd.Series([False] * context.length, dtype="boolean")
        series = context.get_column(config.columns[0])
        return _boolean_series(series).fillna(False)
    if config.builder == "mesh_terms":
        columns = [context.get_column(column) for column in config.columns]
        return _mesh_series(columns)
    if config.builder == "integer":
        if not config.columns:
            return pd.Series([0] * context.length, dtype="Int64")
        series = context.get_column(config.columns[0]).replace("", pd.NA)
        numeric = pd.to_numeric(series, errors="coerce")
        return numeric.astype("Int64").fillna(0)
    if config.builder == "year_column":
        if not config.columns:
            return pd.Series(["" for _ in range(context.length)], dtype="string")
        series = context.get_column(config.columns[0]).replace("", pd.NA)
        numeric = pd.to_numeric(series, errors="coerce")
        return numeric.astype("Int64")
    raise ValueError(f"Unsupported python builder: {config.builder}")


def _build_m_series(context: MPowerQueryContext, config: BuilderConfig) -> pd.Series:
    if config.builder == "column":
        return context.get_column(config.columns[0]) if config.columns else context.get_column("")
    if config.builder == "numeric_identifier":
        result = pd.Series(["" for _ in range(context.length)], dtype="string")
        for column in config.columns:
            candidate = _normalize_numeric_series(context.get_column(column))
            mask = result.eq("") & candidate.ne("")
            result.loc[mask] = candidate.loc[mask]
            if result.ne("").all():
                break
        return result
    if config.builder == "coalesce":
        return _coalesce_columns(context.frame, config.columns)
    if config.builder == "mesh_terms":
        columns = [context.get_column(column) for column in config.columns]
        return _mesh_series(columns)
    if config.builder == "review_flag":
        values = pd.Series([False] * context.length, dtype="boolean")
        for column in config.columns:
            if column not in context.frame.columns:
                continue
            series = _boolean_series(context.frame[column])
            values = values | series.fillna(False)
        return values.fillna(False)
    if config.builder == "publication_year":
        result = pd.Series([pd.NA] * context.length, dtype="Int64")
        for column in config.columns:
            if column not in context.frame.columns:
                continue
            series = context.frame[column]
            if column == "completed":
                candidate = pd.Series([_year_from_date(value) for value in series], dtype="string")
                candidate = candidate.replace("", pd.NA)
                numeric = pd.to_numeric(candidate, errors="coerce")
            else:
                candidate = _string_series(series).replace("", pd.NA)
                numeric = pd.to_numeric(candidate, errors="coerce")
            numeric = numeric.astype("Int64")
            mask = result.isna() & numeric.notna()
            result.loc[mask] = numeric.loc[mask]
            if result.notna().all():
                break
        return result
    if config.builder == "provider_flag":
        provider = config.columns[0] if config.columns else ""
        return context.provider_flag(provider)
    if config.builder == "provider_count":
        providers = config.columns or tuple(PROVIDER_COLUMNS_M.keys())
        flags = [context.provider_flag(provider) for provider in providers]
        if not flags:
            return pd.Series([0] * context.length, dtype="Int64")
        total = pd.DataFrame(flags).T.sum(axis=1)
        return total.astype("Int64")
    raise ValueError(f"Unsupported Power Query builder: {config.builder}")


# ===== Projection builders ===================================================
def _build_projection(frame: pd.DataFrame, crosswalk: Crosswalk, *, side: str) -> pd.DataFrame:
    if side == "python":
        context = PythonDatasetContext(frame)
        builder = _build_python_series
    elif side == "m":
        context = MPowerQueryContext(frame)
        builder = _build_m_series
    else:
        raise ValueError(f"Unknown side {side}")

    result = pd.DataFrame(index=frame.index)
    for field_cfg in crosswalk.fields:
        series = builder(context, getattr(field_cfg, side))
        result[field_cfg.name] = series
    return result


def _canonicalise(df: pd.DataFrame) -> pd.DataFrame:
    canonical = pd.DataFrame(index=df.index)
    for column in df.columns:
        canonical[column] = df[column].map(to_text)
    return canonical


def _build_index(df: pd.DataFrame) -> pd.MultiIndex:
    doc_ids = df["document_chembl_id"].map(to_text)
    pmids = df["primary_pubmed_id"].map(to_text)
    tuples = list(zip(doc_ids, pmids, strict=False))
    return pd.MultiIndex.from_tuples(
        tuples,
        names=("document_chembl_id", "primary_pubmed_id"),
    )


# ===== Metrics ===============================================================
def _group_match_rate(mask: pd.DataFrame, columns: Sequence[str]) -> float:
    if not columns:
        return 1.0
    subset = mask[list(columns)]
    if subset.empty:
        return 1.0
    rowwise = subset.all(axis=1)
    if not len(rowwise):
        return 1.0
    return float(rowwise.sum()) / len(rowwise)


def _per_field_match(mask: pd.DataFrame) -> Mapping[str, float]:
    rates: dict[str, float] = {}
    for column in mask.columns:
        column_mask = mask[column]
        if not len(column_mask):
            rates[column] = 1.0
        else:
            rates[column] = float(column_mask.sum()) / len(column_mask)
    return rates


# ===== Invariants ============================================================
def _invariant_invalid_to_error(
    m_frame: pd.DataFrame,
    python_frame: pd.DataFrame,
) -> Mapping[str, Any]:
    if m_frame.empty:
        return {"passed": True, "violations": 0, "failing_keys": []}

    if "invalid.doi" in m_frame.columns:
        invalid_doi = _boolean_series(m_frame["invalid.doi"]).fillna(False)
    else:
        invalid_doi = pd.Series([False] * len(m_frame), index=m_frame.index, dtype="boolean")

    if "invalid.PMID" in m_frame.columns:
        invalid_pmid = _boolean_series(m_frame["invalid.PMID"]).fillna(False)
    else:
        invalid_pmid = pd.Series([False] * len(m_frame), index=m_frame.index, dtype="boolean")

    if "has_error" in python_frame.columns:
        has_error = _boolean_series(python_frame["has_error"]).fillna(False)
    else:
        has_error = pd.Series([False] * len(python_frame), index=python_frame.index, dtype="boolean")

    if "error_sources" in python_frame.columns:
        error_sources = _string_series(python_frame["error_sources"])
    else:
        error_sources = pd.Series([""] * len(python_frame), index=python_frame.index, dtype="string")

    combined_flag = invalid_doi | invalid_pmid
    violation_mask = combined_flag & ~has_error
    failing_labels: list[Any] = list(m_frame.index[violation_mask])

    for key in combined_flag.index:
        if not bool(combined_flag.loc[key]):
            continue
        tokens = _tokenise_sources(error_sources.loc[key] if key in error_sources.index else "")
        expected_tokens: set[str] = set()
        if bool(invalid_doi.loc[key]):
            expected_tokens |= INVALID_ERROR_TOKEN_MAP["invalid.doi"]
        if bool(invalid_pmid.loc[key]):
            expected_tokens |= INVALID_ERROR_TOKEN_MAP["invalid.PMID"]
        if expected_tokens and tokens.isdisjoint(expected_tokens):
            failing_labels.append(key)

    unique_labels = list(dict.fromkeys(failing_labels))
    sample_index = pd.Index(unique_labels) if unique_labels else pd.Index([])
    sample = _collect_key_sample(sample_index)
    return {
        "passed": not unique_labels,
        "violations": len(unique_labels),
        "failing_keys": sample,
    }


def _invariant_review_merge(
    m_frame: pd.DataFrame,
    python_frame: pd.DataFrame,
) -> Mapping[str, Any]:
    if m_frame.empty:
        return {"passed": True, "violations": 0, "failing_keys": []}

    review = (
        _boolean_series(m_frame["review"]).fillna(False)
        if "review" in m_frame.columns
        else pd.Series([False] * len(m_frame), index=m_frame.index, dtype="boolean")
    )
    doctype = (
        _boolean_series(m_frame["doctype_review"]).fillna(False)
        if "doctype_review" in m_frame.columns
        else pd.Series([False] * len(m_frame), index=m_frame.index, dtype="boolean")
    )
    expected = review | doctype

    if "is_review" in python_frame.columns:
        actual = _boolean_series(python_frame["is_review"]).fillna(False)
    else:
        actual = pd.Series([False] * len(python_frame), index=python_frame.index, dtype="boolean")

    violations = expected & ~actual
    failing_index = pd.Index(m_frame.index[violations])
    sample = _collect_key_sample(failing_index)
    return {
        "passed": int(violations.sum()) == 0,
        "violations": int(violations.sum()),
        "failing_keys": sample,
    }


def _invariant_provider_coverage(
    m_projection: pd.DataFrame,
    python_projection: pd.DataFrame,
    python_raw: pd.DataFrame,
    *,
    provider_columns: Sequence[str],
) -> Mapping[str, Any]:
    if m_projection.empty:
        return {"passed": True, "violations": 0, "failing_keys": []}

    mismatch_flags = pd.Series([False] * len(m_projection), index=m_projection.index, dtype="boolean")

    for provider in provider_columns:
        column = f"has_{provider}"
        if column not in m_projection.columns or column not in python_projection.columns:
            continue
        expected = _boolean_series(m_projection[column]).fillna(False)
        actual = _boolean_series(python_projection[column]).fillna(False)
        mismatch_flags = mismatch_flags | (expected != actual)

    count_expected = pd.to_numeric(
        _string_series(m_projection.get("metadata_source_count", pd.Series([0] * len(m_projection), index=m_projection.index, dtype="string"))).replace("", pd.NA),
        errors="coerce",
    ).fillna(0)
    count_actual = pd.to_numeric(
        _string_series(python_projection.get("metadata_source_count", pd.Series([0] * len(python_projection), index=python_projection.index, dtype="string"))).replace("", pd.NA),
        errors="coerce",
    ).fillna(0)
    count_mismatch = count_expected != count_actual

    coverage_status = _string_series(
        python_raw.get("coverage_status", pd.Series(["" for _ in range(len(python_raw))], index=python_raw.index, dtype="string"))
    ).str.lower()
    has_error = _boolean_series(
        python_raw.get("has_error", pd.Series([False] * len(python_raw), index=python_raw.index, dtype="boolean"))
    ).fillna(False)

    zero_mask = count_expected.eq(0)
    status_violation = zero_mask & ~(coverage_status.isin(["unknown", "failed"]) | (has_error & coverage_status.eq("failed")))

    failing_flags = mismatch_flags | count_mismatch | status_violation
    failing_index = pd.Index(m_projection.index[failing_flags])
    sample = _collect_key_sample(failing_index)
    return {
        "passed": int(failing_flags.sum()) == 0,
        "violations": int(failing_flags.sum()),
        "failing_keys": sample,
    }


# ===== Reporting =============================================================
def _format_percentage(value: float) -> str:
    return f"{value:.4%}"


def _format_invariant(status: Mapping[str, Any]) -> str:
    if not status:
        return "n/a"
    label = "PASS" if status.get("passed") else "FAIL"
    count = status.get("violations", 0)
    return f"{label} ({count} violations)"


def _render_markdown(
    *,
    crosswalk: Crosswalk,
    date_code: str,
    metrics: Mapping[str, Any],
    diff_path: Path | None,
) -> str:
    lines = [f"# QA Report — document postprocessing ({date_code})", ""]
    lines.append(f"- Crosswalk loaded: {crosswalk.version}")
    records = metrics.get("records", {})
    lines.append(
        "- Records: "
        f"python={records.get('python', 0)}, m_like={records.get('m_like', 0)}"
    )
    lines.append("- Match rates:")
    match_rates = metrics.get("match_rates", {})
    for key in ["ids", "preferred", "publication_year", "is_review", "mesh_terms", "has_flags", "coverage"]:
        if key in match_rates:
            lines.append(f"  - {key.replace('_', ' ')}: {_format_percentage(match_rates[key])}")
    invariants = metrics.get("invariants", {})
    lines.append("- Invariants:")
    lines.append(f"  - invalid→has_error: {_format_invariant(invariants.get('invalid_to_error'))}")
    lines.append(f"  - review merge: {_format_invariant(invariants.get('review_merge'))}")
    lines.append(f"  - provider coverage consistent: {_format_invariant(invariants.get('provider_coverage'))}")
    lines.append("- Top 5 mismatch reasons:")
    reasons = metrics.get("top_mismatch_reasons", [])
    if not reasons:
        reasons = ["n/a"]
    for idx, reason in enumerate(reasons[:5], start=1):
        lines.append(f"  {idx}) {reason}")
    lines.append("")
    lines.append(f"Result: {metrics.get('status', 'UNKNOWN')}")
    diff_label = diff_path.name if diff_path else "n/a"
    lines.append(f"Diff: {diff_label}")
    return "\n".join(lines)


# ===== Main QA driver ========================================================
def run_document_postprocessing_check(
    *,
    crosswalk: Crosswalk,
    m_frame: pd.DataFrame,
    python_frame: pd.DataFrame,
    diff_limit: int,
    date_code: str,
    report_dir: Path,
) -> Mapping[str, Any]:
    python_projection = _build_projection(python_frame, crosswalk, side="python")
    m_projection = _build_projection(m_frame, crosswalk, side="m")

    python_canonical = _canonicalise(python_projection)
    m_canonical = _canonicalise(m_projection)

    python_index = _build_index(python_canonical)
    m_index = _build_index(m_canonical)

    python_projection = python_projection.copy()
    m_projection = m_projection.copy()
    python_raw = python_frame.copy()
    m_raw = m_frame.copy()

    python_projection.index = python_index
    m_projection.index = m_index
    python_canonical.index = python_index
    m_canonical.index = m_index
    python_raw.index = python_index
    m_raw.index = m_index

    common_index = python_index.intersection(m_index)
    python_common = python_canonical.reindex(common_index)
    m_common = m_canonical.reindex(common_index)
    python_projection_common = python_projection.reindex(common_index)
    m_projection_common = m_projection.reindex(common_index)
    python_raw_common = python_raw.reindex(common_index)
    m_raw_common = m_raw.reindex(common_index)

    comparison_mask = python_common.eq(m_common)

    per_field = _per_field_match(comparison_mask)
    group_rates: dict[str, float] = {}
    for group, columns in crosswalk.metric_groups().items():
        group_rates[group] = _group_match_rate(comparison_mask, columns)

    python_only = python_index.difference(m_index)
    m_only = m_index.difference(python_index)

    invariants = {
        "invalid_to_error": _invariant_invalid_to_error(
            m_raw_common, python_raw_common
        ),
        "review_merge": _invariant_review_merge(
            m_raw_common, python_raw_common
        ),
        "provider_coverage": _invariant_provider_coverage(
            m_projection_common,
            python_projection_common,
            python_raw_common,
            provider_columns=["chembl", "pubmed", "semantic_scholar", "openalex", "crossref"],
        ),
    }

    mismatched_rows_mask = ~comparison_mask.all(axis=1)
    diff_records: list[dict[str, Any]] = []
    for key in comparison_mask.index[mismatched_rows_mask]:
        row_mask = ~comparison_mask.loc[key]
        differing_columns = row_mask[row_mask].index
        for column in differing_columns:
            diff_records.append(
                {
                    "document_chembl_id": to_text(key[0]),
                    "primary_pubmed_id": to_text(key[1]),
                    "field": column,
                    "python_value": python_common.loc[key, column],
                    "m_value": m_common.loc[key, column],
                }
            )
            if len(diff_records) >= diff_limit:
                break
        if len(diff_records) >= diff_limit:
            break

    diff_path: Path | None = None
    if diff_records:
        diff_df = pd.DataFrame(diff_records)
        report_dir.mkdir(parents=True, exist_ok=True)
        diff_path = report_dir / f"{DIFF_PREFIX}{date_code}{DIFF_SUFFIX}"
        write_csv_deterministic(diff_df, diff_path)

    issues: list[str] = []
    for group, threshold in MATCH_THRESHOLDS.items():
        rate = group_rates.get(group)
        if rate is None:
            continue
        if rate < threshold:
            issues.append(f"Match rate for {group} below threshold ({rate:.4%} < {threshold:.2%})")
    for label, status in invariants.items():
        if not status.get("passed", True):
            issues.append(f"Invariant {label} failed with {status.get('violations', 0)} violations")
    if len(python_only):
        issues.append(f"{len(python_only)} python-only keys detected")
    if len(m_only):
        issues.append(f"{len(m_only)} M-output-only keys detected")

    diff_counts = (~comparison_mask).sum().sort_values(ascending=False)
    top_reasons = [f"{column}: {count}" for column, count in diff_counts.items() if count > 0][:5]

    status_label = "PASS" if not issues else "FAIL"

    metrics = {
        "status": status_label,
        "date_code": date_code,
        "crosswalk_version": crosswalk.version,
        "records": {"python": int(len(python_frame)), "m_like": int(len(m_frame))},
        "missing_rows": {
            "python_only": int(len(python_only)),
            "m_only": int(len(m_only)),
            "python_only_sample": _collect_key_sample(python_only),
            "m_only_sample": _collect_key_sample(m_only),
        },
        "match_rates": group_rates,
        "per_field_match": per_field,
        "invariants": invariants,
        "issues": issues,
        "top_mismatch_reasons": top_reasons,
        "diff_rows": len(diff_records),
    }

    report_dir.mkdir(parents=True, exist_ok=True)
    json_path = report_dir / f"{REPORT_PREFIX}{date_code}{REPORT_JSON_SUFFIX}"
    md_path = report_dir / f"{REPORT_PREFIX}{date_code}{REPORT_MD_SUFFIX}"
    with json_path.open("w", encoding=UTF8_ENCODING) as handle:
        json.dump(metrics, handle, indent=2, ensure_ascii=False)
    markdown = _render_markdown(
        crosswalk=crosswalk,
        date_code=date_code,
        metrics=metrics,
        diff_path=diff_path,
    )
    with md_path.open("w", encoding=UTF8_ENCODING) as handle:
        handle.write(markdown)

    metrics["report_json"] = str(json_path)
    metrics["report_markdown"] = str(md_path)
    metrics["diff_path"] = str(diff_path) if diff_path else None
    return metrics


# ===== CLI ===================================================================
def _build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="QA validation for document post-processing")
    parser.add_argument("--base-path", dest="base_path", default=str(DEFAULT_BASE_PATH))
    parser.add_argument("--out", "--m-output", dest="m_output", required=True)
    parser.add_argument("--preprocessed", dest="python_output", default=None)
    parser.add_argument("--crosswalk", dest="crosswalk_path", default=str(CROSSWALK_PATH))
    parser.add_argument("--diff-limit", dest="diff_limit", type=int, default=DEFAULT_DIFF_LIMIT)
    parser.add_argument("--report-dir", dest="report_dir", default=str(DEFAULT_REPORT_DIR))
    return parser


def main() -> None:
    parser = _build_argument_parser()
    args = parser.parse_args()

    base_path = Path(str(args.base_path).replace("\\", os.sep)).resolve()
    m_path = _resolve_path(base_path, args.m_output)

    if args.python_output:
        python_path = _resolve_path(base_path, args.python_output)
    else:
        python_path = m_path.with_name(f"preprocessed_{m_path.name}")

    crosswalk_path = Path(str(args.crosswalk_path).replace("\\", os.sep))
    if not crosswalk_path.is_absolute():
        crosswalk_path = (Path.cwd() / crosswalk_path).resolve()

    report_dir = Path(str(args.report_dir).replace("\\", os.sep))
    if not report_dir.is_absolute():
        report_dir = (base_path / report_dir).resolve()

    crosswalk = Crosswalk.load(crosswalk_path)
    m_frame = _read_csv(m_path, encoding=CP1252_ENCODING)
    python_frame = _read_csv(python_path, encoding=UTF8_ENCODING)

    date_code = _extract_date_code(python_path)

    metrics = run_document_postprocessing_check(
        crosswalk=crosswalk,
        m_frame=m_frame,
        python_frame=python_frame,
        diff_limit=args.diff_limit,
        date_code=date_code,
        report_dir=report_dir,
    )

    print(json.dumps({"status": metrics.get("status"), "report": metrics.get("report_markdown")}, indent=2))


if __name__ == "__main__":
    main()
