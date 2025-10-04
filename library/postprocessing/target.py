"""Power Query-equivalent target isoform post-processing."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

DEFAULT_ENCODINGS: Sequence[str] = ("utf-8", "utf-8-sig", "cp1252")
SOURCE_COLUMNS: list[str] = [
    "isoform_synonyms",
    "isoform_names",
    "isoform_ids",
    "uniprot_id_primary",
    "target_chembl_id",
]
OUTPUT_COLUMNS: list[str] = ["id", "uniprot_id_primary", "target_chembl_id", "name"]
DEFAULT_INPUT_DIR = Path("data/output")


def _to_text(value: object) -> str:
    """Mirror the Power Query ``ToText`` helper by normalising nulls."""

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if pd.isna(value):  # type: ignore[arg-type]
        return ""
    return str(value)


def _split_pipes(value: object) -> list[str]:
    """Split ``value`` by ``"|"`` removing blanks."""

    text = _to_text(value)
    if not text:
        return []
    parts = [part.strip() for part in text.split("|")]
    return [part for part in parts if part]


def _make_triples(
    names: Iterable[str] | None,
    ids: Iterable[str] | None,
    syns: Iterable[str] | None,
) -> list[dict[str, str | None]]:
    """Align names, ids and synonyms by index padding with ``None``."""

    name_list = list(names or [])
    id_list = list(ids or [])
    syn_list = list(syns or [])
    n = max(len(name_list), len(id_list), len(syn_list))
    if n == 0:
        return []
    triples: list[dict[str, str | None]] = []
    for idx in range(n):
        triples.append(
            {
                "name": name_list[idx] if idx < len(name_list) else None,
                "id": id_list[idx] if idx < len(id_list) else None,
                "synonym": syn_list[idx] if idx < len(syn_list) else None,
            }
        )
    return triples


def _syn_expand(value: object) -> list[str]:
    """Return variants of ``value`` by stripping ``pde`` and ``pld`` substrings."""

    text = _to_text(value).strip().lower()
    if not text:
        return []
    candidates = [text, text.replace("pde", ""), text.replace("pld", "")]
    result: list[str] = []
    seen: set[str] = set()
    for candidate in candidates:
        if candidate and candidate not in seen:
            seen.add(candidate)
            result.append(candidate)
    return result


def _tokenize_synonym(value: object) -> list[str]:
    """Split synonyms on ``":"`` and expand tokens using :func:`_syn_expand`."""

    text = _to_text(value)
    if not text:
        return []
    parts = [segment.strip() for segment in text.split(":")]
    non_empty = [part for part in parts if part]
    variants: list[str] = []
    for part in non_empty:
        for token in _syn_expand(part):
            if token not in variants:
                variants.append(token)
    return variants


def _dedupe_stage1(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"], keep="first"
    )


def _stable_sort(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(
        by=["uniprot_id_primary", "id"], kind="mergesort", na_position="first"
    )


def _dedupe_stage2(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset=["id", "target_chembl_id", "name"], keep="first")


def _dedupe_stage3(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop_duplicates(subset=["id", "name"], keep="first")


def _transform_frame(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    missing = [column for column in SOURCE_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(
            "Input CSV is missing required columns: " + ", ".join(sorted(missing))
        )

    projected = frame[SOURCE_COLUMNS].copy()
    projected["isoform_synonyms"] = projected["isoform_synonyms"].map(
        lambda x: _to_text(x).lower()
    )
    projected["isoform_names"] = projected["isoform_names"].map(
        lambda x: _to_text(x).lower()
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

    rows = projected.explode("triples", ignore_index=True)
    rows = rows[rows["triples"].notna()].reset_index(drop=True)
    if rows.empty:
        expanded = pd.DataFrame(
            columns=[
                "uniprot_id_primary",
                "target_chembl_id",
                "name",
                "id",
                "synonym",
            ]
        )
    else:
        triple_df = pd.DataFrame(rows["triples"].tolist(), columns=["name", "id", "synonym"])
        expanded = pd.concat(
            [
                rows[
                    ["uniprot_id_primary", "target_chembl_id"]
                ].reset_index(drop=True),
                triple_df,
            ],
            axis=1,
        )

    expanded["tokens"] = expanded["synonym"].map(_tokenize_synonym)

    names_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "name"]].copy()
    names_raw["name"] = names_raw["name"].map(lambda x: _to_text(x).strip())
    names_clean = names_raw[~names_raw["name"].isin(["", "n/a", "none"])].copy()

    syns_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "tokens"]].copy()
    syns_exploded = syns_raw.explode("tokens", ignore_index=True)
    if syns_exploded.empty:
        syns_clean = pd.DataFrame(
            columns=["id", "uniprot_id_primary", "target_chembl_id", "name"]
        )
    else:
        syns_exploded["tokens"] = syns_exploded["tokens"].map(
            lambda x: _to_text(x).strip()
        )
        syns_clean = syns_exploded[
            ~syns_exploded["tokens"].isin(["", "n/a", "none"])
        ].copy()
        syns_clean = syns_clean.rename(columns={"tokens": "name"})

    combined = pd.concat([names_clean, syns_clean], ignore_index=True)
    dedup_stage1 = _dedupe_stage1(combined)
    sorted_df = _stable_sort(dedup_stage1)
    dedup_stage2 = _dedupe_stage2(sorted_df)
    result = _dedupe_stage3(dedup_stage2).reset_index(drop=True)

    if result.empty:
        final_df = pd.DataFrame(columns=OUTPUT_COLUMNS)
    else:
        final_df = result[OUTPUT_COLUMNS]

    stats = {
        "combined": len(combined),
        "stage1": len(dedup_stage1),
        "stage2": len(dedup_stage2),
        "stage3": len(final_df),
    }
    return final_df, stats


def _read_csv_with_fallback(
    path: Path, *, encodings: Sequence[str]
) -> tuple[pd.DataFrame, str]:
    last_error: Exception | None = None
    for encoding in encodings:
        try:
            frame = pd.read_csv(path, dtype=object, encoding=encoding)
        except UnicodeDecodeError as exc:
            last_error = exc
            continue
        else:
            return frame, encoding
    if last_error is not None:
        raise last_error
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "unable to decode input")


def _latest_output_file(directory: Path) -> Path:
    candidates = [
        path
        for path in directory.glob("output.target_*.csv")
        if path.is_file()
    ]
    if not candidates:
        raise FileNotFoundError(
            f"No files matching 'output.target_*.csv' found in {directory}"
        )
    return max(candidates, key=lambda item: item.stat().st_mtime)


def process_targets(
    input_csv: str | None = None,
    output_csv: str | None = None,
    verbose: bool = True,
) -> Path:
    """Process target isoform information using the Power Query pipeline."""

    if input_csv is None:
        search_dir = DEFAULT_INPUT_DIR
        input_path = _latest_output_file(search_dir)
    else:
        input_path = Path(input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    if output_csv is None:
        output_path = input_path.with_name(f"isoform.{input_path.name}")
    else:
        output_path = Path(output_csv)

    frame, encoding = _read_csv_with_fallback(input_path, encodings=DEFAULT_ENCODINGS)
    if verbose:
        print(
            f"[target.postprocessing] reading '{input_path}' with encoding {encoding}; rows={len(frame)}"
        )

    transformed, stats = _transform_frame(frame)
    if verbose:
        print(
            f"[target.postprocessing] union rows={stats['combined']}"
        )
        print(
            f"[target.postprocessing] dedup stage1 {stats['combined']} -> {stats['stage1']}"
        )
        print(
            f"[target.postprocessing] dedup stage2 {stats['stage1']} -> {stats['stage2']}"
        )
        print(
            f"[target.postprocessing] dedup stage3 {stats['stage2']} -> {stats['stage3']}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    transformed.to_csv(output_path, index=False, encoding="utf-8")

    if verbose:
        unique_ids = transformed["id"].nunique(dropna=True) if not transformed.empty else 0
        print(
            f"[target.postprocessing] wrote '{output_path}' with {stats['stage3']} rows; unique ids={unique_ids}"
        )
    return output_path
