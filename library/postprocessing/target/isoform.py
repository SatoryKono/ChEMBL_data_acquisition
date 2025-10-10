"""Isoform post-processing helpers for target exports.

This module mirrors the Power Query workbook used for the isoform expansion
stage of the targets pipeline. The implementation is intentionally verbose to
match the original M transformations byte-for-byte and to keep the resulting
CSV deterministic for regression testing.
"""

from __future__ import annotations

import re
import sys
import warnings
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from ..helpers import normalise_export_basename


def _empty_like(index: pd.Index) -> pd.Series:
    """Return an empty object-typed series aligned with *index*."""

    if len(index):
        return pd.Series([""] * len(index), index=index, dtype=object)
    return pd.Series(dtype=object)


def _project_source_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` projected onto :data:`_SOURCE_COLUMNS` with fallbacks."""

    projected = frame.reindex(columns=list(_SOURCE_COLUMNS)).copy()
    for column in _SOURCE_COLUMNS:
        if column in frame.columns:
            projected[column] = frame[column].astype(object)
            continue

        fallback_series: pd.Series | None = None
        for candidate in _SOURCE_FALLBACKS.get(
            column, ()
        ):  # pragma: no branch - small tuple
            if candidate in frame.columns:
                fallback_series = frame[candidate].astype(object)
                break
            if candidate in projected.columns and candidate != column:
                existing = projected[candidate]
                if not existing.replace("", pd.NA).isna().all():
                    fallback_series = existing.astype(object)
                    break

        if fallback_series is None:
            projected[column] = _empty_like(projected.index)
        else:
            projected[column] = fallback_series.reindex(projected.index).astype(object)

    return projected


# Ordered list of encodings attempted when reading the aggregated targets CSV.
_DEFAULT_ENCODINGS: tuple[str, ...] = ("utf-8", "utf-8-sig", "cp1252")

# Directory searched when ``process_targets`` is invoked without an explicit
# path. The Power Query workflow consumed the canonical exports emitted under
# ``data/output`` so we keep the same convention here.
_DEFAULT_SEARCH_DIR = Path("data/output")

# Accepted filename patterns for aggregated target exports.
_INPUT_NAME_RULES: tuple[tuple[re.Pattern[str], str], ...] = (
    (
        re.compile(r"output\.target_\d{8}(?:_[A-Za-z0-9]+)*(?:_normalized)?\.csv\Z"),
        "output.target_<YYYYMMDD>[<_suffixes>][_normalized].csv",
    ),
    (
        re.compile(r"output\.targets_\d{8}(?:_[A-Za-z0-9]+)*(?:_normalized)?\.csv\Z"),
        "output.targets_<YYYYMMDD>[<_suffixes>][_normalized].csv",
    ),
    (
        re.compile(r"output\.targets(?:_[A-Za-z0-9]+)*(?:_normalized)?\.csv\Z"),
        "output.targets[_<suffixes>][_normalized].csv",
    ),
    (
        re.compile(r"targets_\d{8}(?:_[A-Za-z0-9]+)*(?:_normalized)?\.csv\Z"),
        "targets_<YYYYMMDD>[<_suffixes>][_normalized].csv",
    ),
    (
        re.compile(r"out_uniprot(?:_[A-Za-z0-9]+)*\.csv\Z"),
        "out_uniprot[_<suffixes>].csv",
    ),
    (
        re.compile(r"out(?:_[A-Za-z0-9]+)*\.csv\Z"),
        "out[_<suffixes>].csv",
    ),
)


def _matches_expected_input_name(filename: str) -> bool:
    """Return ``True`` when ``filename`` matches a supported input pattern."""

    return any(pattern.match(filename) for pattern, _ in _INPUT_NAME_RULES)


def _supported_patterns_text() -> str:
    """Return a human-readable list of accepted filename templates."""

    return ", ".join(description for _, description in _INPUT_NAME_RULES)


def _current_default_search_dir() -> Path:
    """Return the current default search directory with package overrides."""

    package = sys.modules.get("library.postprocessing.target")
    if package is not None and hasattr(package, "_DEFAULT_SEARCH_DIR"):
        override = package._DEFAULT_SEARCH_DIR
        if override is not None:
            return Path(override)
    return _DEFAULT_SEARCH_DIR


# Columns projected from the aggregated target table before isoform expansion.
_SOURCE_COLUMNS: tuple[str, ...] = (
    "isoform_synonyms",
    "isoform_names",
    "isoform_ids",
    "uniprot_id_primary",
    "target_chembl_id",
)

# Alternative column names consulted when projecting the input frame.  The
# UniProt-only export, for instance, lacks ``target_chembl_id`` and
# ``uniprot_id_primary`` because the identifiers are represented as
# ``uniProtkbId`` (and friends).  The fallbacks ensure that these shapes can be
# processed without raising ``KeyError`` while keeping legacy behaviour for the
# canonical aggregated target table.
_SOURCE_FALLBACKS: dict[str, tuple[str, ...]] = {
    "uniprot_id_primary": (
        "uniprot_id",
        "uniprotkb_Id",
        "uniprotkb_id",
        "uniProtkbId",
        "uniProtkbIdFallback",
    ),
    "target_chembl_id": (
        "target_id",
        "chembl_id",
        "target",
        "uniprot_id_primary",
        "uniprot_id",
        "uniprotkb_Id",
        "uniProtkbId",
        "uniProtkbIdFallback",
    ),
}

_SOURCE_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "uniprot_id_primary": (
        "uniprot_id",
        "primary_accession",
        "uniprot_accession",
        "accession",
    ),
    "target_chembl_id": (
        "chembl_id",
        "target_id",
        "target_chemblid",
    ),
}

# Column order of the emitted CSV artefact.
_OUTPUT_COLUMNS: tuple[str, ...] = (
    "id",
    "uniprot_id_primary",
    "target_chembl_id",
    "name",
)


def _normalize_column_name(column: str) -> str:
    """Return a normalised representation of ``column`` for fuzzy matching."""

    return re.sub(r"[^a-z0-9]", "", _to_text(column).lower())


def _candidate_source_names(canonical: str) -> tuple[str, ...]:
    """Return the ordered list of accepted names for ``canonical``."""

    seen: list[str] = []
    for name in (
        (canonical,)
        + _SOURCE_COLUMN_ALIASES.get(canonical, ())
        + _SOURCE_FALLBACKS.get(canonical, ())
    ):
        if name not in seen:
            seen.append(name)
    return tuple(seen)


def _format_alias_hint(canonical: str) -> str:
    """Return a human-readable description of accepted names."""

    candidates = _candidate_source_names(canonical)
    if len(candidates) == 1:
        return canonical
    aliases = ", ".join(candidates[1:])
    return f"{canonical} (aliases: {aliases})"


def _resolve_source_columns(frame: pd.DataFrame) -> dict[str, str]:
    """Resolve the source column names accounting for legacy aliases."""

    normalized_lookup: dict[str, list[str]] = {}
    for column in frame.columns:
        normalized_lookup.setdefault(_normalize_column_name(column), []).append(column)

    resolved: dict[str, str] = {}
    missing: list[str] = []

    for canonical in _SOURCE_COLUMNS:
        match: str | None = None
        candidates = _candidate_source_names(canonical)

        for candidate in candidates:  # pragma: no branch - small tuple
            if candidate in frame.columns:
                match = candidate
                break

        if match is None:
            candidate_norms = {_normalize_column_name(name) for name in candidates}
            for norm in candidate_norms:
                candidates_from_norm = normalized_lookup.get(norm)
                if candidates_from_norm:
                    match = sorted(candidates_from_norm)[0]
                    break

        if match is None:
            missing.append(canonical)
        else:
            resolved[canonical] = match

    if missing:
        alias_hints = ", ".join(_format_alias_hint(name) for name in missing)
        available = ", ".join(frame.columns)
        raise KeyError(
            "Required columns missing from isoform export: "
            + f"{missing}. "
            + f"Available columns: {available}. "
            + f"Accepted names: {alias_hints}"
        )

    return resolved


def _to_text(value: Any) -> str:
    """Replicate the Power Query ``ToText`` helper."""

    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", errors="ignore")
    if isinstance(value, (int, float)) and pd.isna(value):
        return ""
    return str(value)


def _split_pipes(value: Any) -> list[str]:
    """Split values by ``"|"`` trimming whitespace and dropping blanks."""

    text = _to_text(value)
    if not text:
        return []
    parts = [segment.strip() for segment in text.split("|")]
    return [segment for segment in parts if segment]


def _make_triples(
    names: Sequence[str], ids: Sequence[str], synonyms: Sequence[str]
) -> list[dict[str, str | None]]:
    """Align isoform names/ids/synonyms into indexed triplets."""

    n = max(len(names), len(ids), len(synonyms))
    if n == 0:
        return []
    triples: list[dict[str, str | None]] = []
    for idx in range(n):
        triples.append(
            {
                "name": names[idx] if idx < len(names) else None,
                "id": ids[idx] if idx < len(ids) else None,
                "synonym": synonyms[idx] if idx < len(synonyms) else None,
            }
        )
    return triples


def _syn_expand(value: Any) -> list[str]:
    """Expand a synonym token according to the Power Query rules."""

    token = _to_text(value).lower().strip()
    if not token:
        return []
    candidates: list[str] = []
    for variant in (token, token.replace("pde", ""), token.replace("pld", "")):
        variant = variant.strip()
        if variant and variant not in candidates:
            candidates.append(variant)
    return candidates


def _tokenize_synonym(value: Any) -> list[str]:
    """Tokenise isoform synonyms by ``":"`` with ``SynExpand`` variants."""

    text = _to_text(value)
    if not text:
        return []
    parts = [segment.strip() for segment in text.split(":")]
    tokens: list[str] = []
    for part in parts:
        if not part:
            continue
        for variant in _syn_expand(part):
            if variant not in tokens:
                tokens.append(variant)
    return tokens


@dataclass(frozen=True)
class _TransformationResult:
    result: pd.DataFrame
    combined: pd.DataFrame
    dedup_stage1: pd.DataFrame
    sorted_stage: pd.DataFrame
    dedup_stage2: pd.DataFrame


def _transform(frame: pd.DataFrame) -> _TransformationResult:
    """Apply the isoform expansion stages mirroring the Power Query steps."""

    def _empty_result() -> _TransformationResult:
        empty = pd.DataFrame(columns=list(_OUTPUT_COLUMNS))
        return _TransformationResult(
            result=empty.copy(),
            combined=empty.copy(),
            dedup_stage1=empty.copy(),
            sorted_stage=empty.copy(),
            dedup_stage2=empty.copy(),
        )

    if frame.empty:
        return _empty_result()

    try:
        resolved = _resolve_source_columns(frame)
    except KeyError as exc:
        warnings.warn(
            (
                "Isoform post-processing could not resolve required columns; "
                "emitting an empty result instead. Details: %s"
            )
            % exc,
            UserWarning,
            stacklevel=2,
        )
        return _empty_result()
    aligned = frame.copy()
    for canonical, source in resolved.items():
        if source in aligned.columns:
            aligned[canonical] = aligned[source].astype(object)
        else:  # pragma: no cover - defensive guard for mutated inputs
            aligned[canonical] = _empty_like(aligned.index)

    projected = _project_source_columns(aligned)
    if projected.empty:
        return _empty_result()
    projected = projected.loc[:, list(_SOURCE_COLUMNS)].copy()

    projected["isoform_synonyms"] = projected["isoform_synonyms"].map(
        lambda value: _to_text(value).lower()
    )
    projected["isoform_names"] = projected["isoform_names"].map(
        lambda value: _to_text(value).lower()
    )

    projected["isoform_synonyms"] = projected["isoform_synonyms"].map(_split_pipes)
    projected["isoform_names"] = projected["isoform_names"].map(_split_pipes)
    projected["isoform_ids"] = projected["isoform_ids"].map(_split_pipes)

    projected["triples"] = projected.apply(
        lambda row: _make_triples(
            row["isoform_names"], row["isoform_ids"], row["isoform_synonyms"]
        ),
        axis=1,
    )

    expanded = projected.explode("triples", ignore_index=True)
    expanded = expanded[expanded["triples"].notna()].copy()
    triple_cols = expanded["triples"].apply(pd.Series)
    expanded = pd.concat([expanded.drop(columns=["triples"]), triple_cols], axis=1)

    expanded["tokens"] = expanded["synonym"].map(_tokenize_synonym)

    names_raw = expanded[
        ["id", "uniprot_id_primary", "target_chembl_id", "name"]
    ].copy()
    names_raw["name"] = names_raw["name"].map(lambda value: _to_text(value).strip())
    names_clean = names_raw[~names_raw["name"].isin({"", "n/a", "none"})].reset_index(
        drop=True
    )

    syns_raw = expanded[
        ["id", "uniprot_id_primary", "target_chembl_id", "tokens"]
    ].copy()
    syns_exploded = syns_raw.explode("tokens", ignore_index=True)
    syns_exploded["tokens"] = syns_exploded["tokens"].map(
        lambda value: _to_text(value).strip()
    )
    syns_clean = syns_exploded[
        ~syns_exploded["tokens"].isin({"", "n/a", "none"})
    ].rename(columns={"tokens": "name"})

    combined = pd.concat([names_clean, syns_clean], ignore_index=True)

    dedup_stage1 = combined.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"]
    )

    sorted_stage = dedup_stage1.sort_values(
        by=["uniprot_id_primary", "id"], kind="mergesort", na_position="first"
    ).reset_index(drop=True)

    dedup_stage2 = sorted_stage.drop_duplicates(
        subset=["id", "target_chembl_id", "name"], keep="first"
    )

    result = dedup_stage2.drop_duplicates(subset=["id", "name"], keep="first")
    result = result.loc[:, list(_OUTPUT_COLUMNS)].reset_index(drop=True)

    return _TransformationResult(
        result=result,
        combined=combined,
        dedup_stage1=dedup_stage1,
        sorted_stage=sorted_stage,
        dedup_stage2=dedup_stage2,
    )


def _read_csv(path: Path, *, encodings: Iterable[str]) -> tuple[pd.DataFrame, str]:
    """Read ``path`` trying encodings sequentially until success."""

    last_error: Exception | None = None
    for encoding in encodings:
        try:
            frame = pd.read_csv(
                path,
                dtype=str,
                encoding=encoding,
                sep=",",
                keep_default_na=False,
            )
        except UnicodeDecodeError as exc:
            last_error = exc
            continue
        return frame, encoding
    if last_error is not None:  # pragma: no cover - defensive guard
        raise last_error
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "unable to decode input")


def _resolve_input_path(input_csv: str | Path | None) -> Path:
    """Resolve the source CSV path following the documented conventions."""

    if input_csv is not None:
        candidate = Path(input_csv)
        if candidate.is_dir():
            matches = [
                path
                for path in candidate.glob("*.csv")
                if _matches_expected_input_name(path.name)
            ]
            if not matches:
                raise FileNotFoundError(
                    "No supported target exports found under "
                    f"{candidate} (expected patterns: {_supported_patterns_text()})"
                )
            return max(matches, key=lambda path: (path.stat().st_mtime, path.name))
        if not candidate.exists():
            raise FileNotFoundError(candidate)
        if not _matches_expected_input_name(candidate.name):
            if candidate.suffix.lower() != ".csv":
                warnings.warn(
                    (
                        "Input file '%s' does not use the canonical '.csv' extension "
                        "(%s); attempting to proceed because an explicit path was "
                        "provided."
                    )
                    % (candidate.name, _supported_patterns_text()),
                    UserWarning,
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    (
                        "Input file '%s' does not match the canonical target export "
                        "naming conventions (%s); proceeding because an explicit path "
                        "was provided."
                    )
                    % (candidate.name, _supported_patterns_text()),
                    UserWarning,
                    stacklevel=2,
                )
        return candidate

    search_dir = _current_default_search_dir()
    if not search_dir.exists():
        raise FileNotFoundError(
            "No input CSV supplied and default search directory does not exist"
        )
    matches = [
        path
        for path in search_dir.glob("*.csv")
        if _matches_expected_input_name(path.name)
    ]
    if not matches:
        raise FileNotFoundError(
            "No supported target exports found under "
            f"{search_dir} (expected patterns: {_supported_patterns_text()})"
        )
    return max(matches, key=lambda path: (path.stat().st_mtime, path.name))


def _resolve_output_path(input_path: Path, output_csv: str | None) -> Path:
    if output_csv is not None:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return output_path

    base_name = normalise_export_basename(input_path)
    return input_path.with_name(f"isoform.{base_name}")


def _stringify_for_csv(value: Any) -> str:
    text = _to_text(value)
    return text.strip()


def process_targets(
    input_csv: str | None = None,
    output_csv: str | None = None,
    verbose: bool = True,
) -> Path:
    """Run the isoform post-processing pipeline on canonical target exports."""

    input_path = _resolve_input_path(input_csv)

    output_path = _resolve_output_path(input_path, output_csv)

    frame, encoding_used = _read_csv(input_path, encodings=_DEFAULT_ENCODINGS)
    if verbose:
        print(f"[isoform] read {input_path} using encoding {encoding_used}")

    transform = _transform(frame)

    if verbose:
        print(f"[isoform] combined rows before dedup: {len(transform.combined)}")
        print("[isoform] rows after stage1: " f"{len(transform.dedup_stage1)}")
        print("[isoform] rows after stage2: " f"{len(transform.dedup_stage2)}")
        print("[isoform] final rows: " f"{len(transform.result)}")

    prepared = transform.result.copy()
    for column in _OUTPUT_COLUMNS:
        prepared[column] = prepared[column].map(_stringify_for_csv)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(
        output_path,
        index=False,
        encoding="utf-8",
        lineterminator="\n",
    )
    if verbose:
        print(f"[isoform] wrote {output_path}")
    return output_path


__all__ = [
    "process_targets",
    "_split_pipes",
    "_make_triples",
    "_syn_expand",
    "_tokenize_synonym",
    "_transform",
    "_matches_expected_input_name",
    "_DEFAULT_SEARCH_DIR",
    "_DEFAULT_ENCODINGS",
    "_INPUT_NAME_RULES",
    "_SOURCE_COLUMNS",
    "_OUTPUT_COLUMNS",
]
