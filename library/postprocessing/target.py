"""Isoform post-processing helpers for target exports.

This module mirrors the Power Query workbook used for the isoform expansion
stage of the targets pipeline. The implementation is intentionally verbose to
match the original M transformations byte-for-byte and to keep the resulting
CSV deterministic for regression testing.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

# Ordered list of encodings attempted when reading the aggregated targets CSV.
_DEFAULT_ENCODINGS: tuple[str, ...] = ("utf-8", "utf-8-sig", "cp1252")

# Directory searched when ``process_targets`` is invoked without an explicit
# path. The Power Query workflow consumed the canonical exports emitted under
# ``data/output`` so we keep the same convention here.
_DEFAULT_SEARCH_DIR = Path("data/output")

# Columns projected from the aggregated target table before isoform expansion.
_SOURCE_COLUMNS: tuple[str, ...] = (
    "isoform_synonyms",
    "isoform_names",
    "isoform_ids",
    "uniprot_id_primary",
    "target_chembl_id",
)

# Column order of the emitted CSV artefact.
_OUTPUT_COLUMNS: tuple[str, ...] = (
    "id",
    "uniprot_id_primary",
    "target_chembl_id",
    "name",
)


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

    projected = frame.loc[:, list(_SOURCE_COLUMNS)].copy()

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

    names_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "name"]].copy()
    names_raw["name"] = names_raw["name"].map(lambda value: _to_text(value).strip())
    names_clean = names_raw[
        ~names_raw["name"].isin({"", "n/a", "none"})
    ].reset_index(drop=True)

    syns_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "tokens"]].copy()
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
            matches = list(candidate.glob("output.target_*.csv"))
            if not matches:
                raise FileNotFoundError(
                    f"No output.target_*.csv files found under {candidate}"
                )
            return max(matches, key=lambda path: path.stat().st_mtime)
        if not candidate.exists():
            raise FileNotFoundError(candidate)
        return candidate

    search_dir = _DEFAULT_SEARCH_DIR
    if not search_dir.exists():
        raise FileNotFoundError(
            "No input CSV supplied and default search directory does not exist"
        )
    matches = list(search_dir.glob("output.target_*.csv"))
    if not matches:
        raise FileNotFoundError(
            f"No output.target_*.csv files found under {search_dir}" 
        )
    return max(matches, key=lambda path: path.stat().st_mtime)


def _resolve_output_path(input_path: Path, output_csv: str | None) -> Path:
    if output_csv is not None:
        output_path = Path(output_csv)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        return output_path
    return input_path.with_name(f"isoform.{input_path.name}")


def _stringify_for_csv(value: Any) -> str:
    text = _to_text(value)
    return text.strip()


def process_targets(
    input_csv: str | None = None,
    output_csv: str | None = None,
    verbose: bool = True,
) -> Path:
    """Run the isoform post-processing pipeline on ``output.target_*.csv``."""

    input_path = _resolve_input_path(input_csv)
    if not input_path.name.startswith("output.target_"):
        raise ValueError("Input file must match pattern output.target_*.csv")

    output_path = _resolve_output_path(input_path, output_csv)

    frame, encoding_used = _read_csv(input_path, encodings=_DEFAULT_ENCODINGS)
    if verbose:
        print(f"[isoform] read {input_path} using encoding {encoding_used}")

    transform = _transform(frame)

    if verbose:
        print(f"[isoform] combined rows before dedup: {len(transform.combined)}")
        print(
            "[isoform] rows after stage1: "
            f"{len(transform.dedup_stage1)}"
        )
        print(
            "[isoform] rows after stage2: "
            f"{len(transform.dedup_stage2)}"
        )
        print(
            "[isoform] final rows: "
            f"{len(transform.result)}"
        )

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
]
