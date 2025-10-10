"""IUPHAR post-processing helpers reproducing the legacy Power Query workbook."""

from __future__ import annotations

import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from library.common.csv_utils import write_csv_deterministic
from library.common.log import logger

from .helpers import normalise_export_basename, normalise_text, read_csv_with_fallbacks

__all__ = ["process_iuphar_targets", "IUPHARPostProcessingError"]

_DEFAULT_SEARCH_DIR = Path("data/output")
_TARGET_PATTERN = re.compile(r"output\.target_\d{8}\.csv\Z")


class IUPHARPostProcessingError(RuntimeError):
    """Raised when the IUPHAR helper cannot proceed."""


@dataclass(frozen=True)
class SynonymStatistics:
    before: int = 0
    after: int = 0


_DROP_COLUMNS: tuple[str, ...] = (
    "component_synonym_ids",
    "component_type_raw",
    "component_sequence",
    "component_structures",
)

_RENAME_COLUMNS: dict[str, str] = {
    "GuidetoPHARMACOLOGY": "guidetopharmacology_id",
}

_REQUIRED_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "GuidetoPHARMACOLOGY",
    "iuphar_target_id",
    "iuphar_family_id",
    "iuphar_type",
    "iuphar_class",
    "iuphar_subclass",
    "iuphar_chain",
    "iuphar_name",
    "gtop_synonyms",
)

_OPTIONAL_COLUMNS: tuple[str, ...] = ("component_description",)

_OUTPUT_COLUMNS: tuple[str, ...] = (
    "target_chembl_id",
    "guidetopharmacology_id",
    "iuphar_target_id",
    "iuphar_family_id",
    "iuphar_type",
    "iuphar_class",
    "iuphar_subclass",
    "iuphar_chain",
    "iuphar_name",
    "iuphar_synonyms",
)


def _current_default_search_dir() -> Path:
    package = sys.modules.get(__name__)
    if package is not None and hasattr(package, "_DEFAULT_SEARCH_DIR"):
        override = package._DEFAULT_SEARCH_DIR
        if override is not None:
            return Path(override)
    return _DEFAULT_SEARCH_DIR


def _matches_target_filename(name: str) -> bool:
    return bool(_TARGET_PATTERN.match(name))


def _latest_target_file(search_dir: Path) -> Path:
    candidates = sorted(
        (
            path
            for path in search_dir.iterdir()
            if path.is_file() and _matches_target_filename(path.name)
        ),
        key=lambda item: item.name,
    )
    if not candidates:
        raise IUPHARPostProcessingError(
            f"No target exports matching 'output.target_YYYYMMDD.csv' found in {search_dir!s}"
        )
    return candidates[-1]


def _normalise_column_name(name: str) -> str:
    """Return ``name`` stripped of whitespace and UTF-8 BOM markers."""

    if isinstance(name, str):
        # ``str.strip`` alone does not remove the BOM prefix, therefore we
        # explicitly lstrip it before trimming whitespace.  The helper keeps the
        # original casing because downstream code relies on exact header names.
        return name.lstrip("\ufeff").strip()
    return name


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return ``df`` with consistently normalised column headers."""

    if not df.columns.empty:
        df = df.rename(columns=_normalise_column_name)
    return df


def _ensure_required_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, tuple[str, ...]]:
    missing = tuple(column for column in _REQUIRED_COLUMNS if column not in df.columns)
    if not missing:
        return df, ()

    logger.warning(
        "iuphar_postprocess_missing_columns",
        missing=list(missing),
        available=list(df.columns),
    )

    formatted = ", ".join(missing)
    raise IUPHARPostProcessingError(
        f"Input CSV is missing required columns: {formatted}"
    )


def _ensure_optional_columns(df: pd.DataFrame) -> pd.DataFrame:
    missing = [column for column in _OPTIONAL_COLUMNS if column not in df.columns]
    if not missing:
        return df
    df = df.copy()
    for column in missing:
        df[column] = ""
    return df


def _clean_brackets(value: str) -> str:
    text = re.sub(r"\[[^\]]*\]", "", value)
    text = re.sub(r"\([^\)]*\)", "", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _split_pipe(value: str) -> list[str]:
    tokens = []
    for token in value.split("|"):
        token = token.strip()
        if token:
            tokens.append(token)
    return tokens


def _parse_component_descriptions(value: str) -> list[str]:
    text = value.strip()
    if not text or text in {"-", "null", "NULL"}:
        return []
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        tokens = re.split(r"[|;]", text)
        return [item.strip() for item in tokens if item.strip()]
    records = payload if isinstance(payload, list) else [payload]
    result: list[str] = []
    for record in records:
        if isinstance(record, dict):
            for key in ("component_description", "description", "name"):
                raw = record.get(key)
                if isinstance(raw, str) and raw.strip():
                    result.append(raw.strip())
        elif isinstance(record, str) and record.strip():
            result.append(record.strip())
    return result


def _collect_synonym_tokens(row: pd.Series[Any]) -> tuple[list[str], list[str]]:
    raw_tokens: list[str] = []
    sources: tuple[str | None, ...] = (
        (
            normalise_text(row.get("gtop_synonyms"))
            if "gtop_synonyms" in row and pd.notnull(row["gtop_synonyms"])
            else None
        ),
        (
            normalise_text(row.get("synonyms"))
            if "synonyms" in row and pd.notnull(row["synonyms"])
            else None
        ),
    )

    for source in sources:
        if source:
            raw_tokens.extend(_split_pipe(source))
    component_desc = normalise_text(row.get("component_description"))
    if component_desc:
        raw_tokens.extend(_parse_component_descriptions(component_desc))
    name_value = normalise_text(row.get("iuphar_name"))
    if name_value:
        raw_tokens.append(name_value)
    cleaned_tokens: list[str] = []
    for token in raw_tokens:
        text = normalise_text(token)
        if not text:
            continue
        text = _clean_brackets(text)
        text = text.lower()
        if text:
            cleaned_tokens.append(text)
    deduped: list[str] = []
    seen: set[str] = set()
    for token in cleaned_tokens:
        if token not in seen:
            seen.add(token)
            deduped.append(token)
    return cleaned_tokens, deduped


def _apply_synonym_processing(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, SynonymStatistics]:
    before = 0
    after = 0
    deduped_rows: list[list[str]] = []
    for _, row in df.iterrows():
        cleaned, deduped = _collect_synonym_tokens(row)
        before += len(cleaned)
        after += len(deduped)
        deduped_rows.append(deduped)
    df = df.copy()
    df["iuphar_synonyms"] = ["|".join(tokens) for tokens in deduped_rows]
    return df, SynonymStatistics(before=before, after=after)


def _drop_helper_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    drop_candidates = [column for column in _DROP_COLUMNS if column in df.columns]
    if not drop_candidates:
        return df, 0
    return df.drop(columns=drop_candidates, errors="ignore"), len(drop_candidates)


def _prepare_output(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {old: new for old, new in _RENAME_COLUMNS.items() if old in df.columns}
    df = df.rename(columns=rename_map)
    for column in _OUTPUT_COLUMNS:
        if column not in df.columns:
            df[column] = ""
    return df.loc[:, list(_OUTPUT_COLUMNS)]


def process_iuphar_targets(
    input_csv: str | Path | None = None,
    *,
    output_csv: str | Path | None = None,
    verbose: bool = False,
) -> Path:
    """Process the canonical target export to reproduce the legacy IUPHAR helper."""

    if input_csv is None:
        search_dir = _current_default_search_dir()
        input_path = _latest_target_file(Path(search_dir))
    else:
        input_path = Path(input_csv)
    if output_csv is None:
        base = normalise_export_basename(input_path)
        output_path = input_path.with_name(f"IUPHAR.{base}").with_suffix(".csv")
    else:
        output_path = Path(output_csv)

    df = read_csv_with_fallbacks(input_path)
    df = _normalise_columns(df)
    df, missing_required = _ensure_required_columns(df)
    df = _ensure_optional_columns(df)

    input_rows = len(df)
    df, dropped_columns = _drop_helper_columns(df)
    df, stats = _apply_synonym_processing(df)
    df = _prepare_output(df)
    output_rows = len(df)

    if missing_required:
        df = df.iloc[0:0]
        if verbose:
            logger.info(
                "iuphar_postprocess: emitted empty result due to missing columns",
                extra={
                    "missing_columns": list(missing_required),
                    "input_path": str(input_path),
                },
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Write with stable column order and LF newlines without BOM to match snapshots
    write_csv_deterministic(
        df,
        output_path,
        col_order=_OUTPUT_COLUMNS,
        key_cols=("target_chembl_id",),
        encoding="utf-8",
    )

    if verbose:
        logger.info("iuphar_postprocess: input=%s output=%s", input_path, output_path)
        logger.info(
            "iuphar_postprocess: input_rows=%d output_rows=%d dropped_columns=%d",
            input_rows,
            output_rows,
            dropped_columns,
        )
        logger.info(
            "iuphar_postprocess: synonym_tokens_before=%d synonym_tokens_after=%d",
            stats.before,
            stats.after,
        )

    return output_path
