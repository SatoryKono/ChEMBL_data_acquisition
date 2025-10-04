"""Post-processing helpers for target isoform tables.

The transformations mirror the Power Query (M) logic historically used to
derive isoform name/synonym triples from the primary target export. The goal is
to provide a deterministic, bit-for-bit equivalent implementation in Python.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import pandas as pd

__all__ = [
    "process_targets",
    "_split_pipes",
    "_make_triples",
    "_syn_expand",
    "_tokenize_synonym",
    "_deduplicate_stage1",
    "_sort_for_stage2",
    "_deduplicate_stage2",
    "_deduplicate_final",
]


# ===== Parameters ===========================================================

TARGET_INPUT_COLUMNS: Sequence[str] = (
    "isoform_synonyms",
    "isoform_names",
    "isoform_ids",
    "uniprot_id_primary",
    "target_chembl_id",
)

FINAL_COLUMN_ORDER: Sequence[str] = (
    "id",
    "uniprot_id_primary",
    "target_chembl_id",
    "name",
)

CSV_ENCODING_PREFERENCE: Sequence[str] = ("utf-8", "utf-8-sig", "cp1252")

DEFAULT_SEARCH_DIR = Path("output")


# ===== Helper utilities =====================================================


def _to_text(value: object) -> str:
    """Replicate Power Query's ``ToText`` helper."""

    if value is None:
        return ""
    if isinstance(value, float) and pd.isna(value):  # type: ignore[arg-type]
        return ""
    if pd.isna(value):  # type: ignore[arg-type]
        return ""
    return str(value)


def _split_pipes(value: object) -> list[str]:
    """Split pipe-delimited strings removing whitespace and empty entries."""

    text = _to_text(value)
    if not text:
        return []
    parts = [segment.strip() for segment in text.split("|")]
    return [part for part in parts if part]


def _make_triples(
    names: Sequence[str], ids: Sequence[str], synonyms: Sequence[str]
) -> list[dict[str, str | None]]:
    """Align name/id/synonym lists by index padding with ``None``.

    Parameters
    ----------
    names, ids, synonyms:
        Lists obtained after splitting the respective isoform columns.
    """

    n_names = len(names)
    n_ids = len(ids)
    n_synonyms = len(synonyms)
    size = max(n_names, n_ids, n_synonyms)
    if size == 0:
        return []
    triples: list[dict[str, str | None]] = []
    for index in range(size):
        triples.append(
            {
                "name": names[index] if index < n_names else None,
                "id": ids[index] if index < n_ids else None,
                "synonym": synonyms[index] if index < n_synonyms else None,
            }
        )
    return triples


def _syn_expand(token: object) -> list[str]:
    """Return token variants with ``pde``/``pld`` removed when present."""

    base = _to_text(token).strip().lower()
    if not base:
        return []
    variants = [base, base.replace("pde", ""), base.replace("pld", "")]
    unique: list[str] = []
    seen: set[str] = set()
    for variant in variants:
        if not variant:
            continue
        if variant in seen:
            continue
        unique.append(variant)
        seen.add(variant)
    return unique


def _tokenize_synonym(value: object) -> list[str]:
    """Split synonyms on ``:`` and expand ``pde``/``pld`` variants."""

    raw = _to_text(value)
    if not raw:
        return []
    parts = [segment.strip() for segment in raw.split(":")]
    non_empty = [segment for segment in parts if segment]
    tokens: list[str] = []
    for segment in non_empty:
        for variant in _syn_expand(segment):
            if variant not in tokens:
                tokens.append(variant)
    return tokens


def _deduplicate_stage1(df: pd.DataFrame) -> pd.DataFrame:
    """First deduplication stage on four columns."""

    subset = ["id", "name", "target_chembl_id", "uniprot_id_primary"]
    return df.drop_duplicates(subset=subset, keep="first").reset_index(drop=True)


def _sort_for_stage2(df: pd.DataFrame) -> pd.DataFrame:
    """Deterministic sort prior to the second deduplication stage."""

    if df.empty:
        return df.copy()
    return df.sort_values(
        by=["uniprot_id_primary", "id"],
        kind="mergesort",
        na_position="first",
    ).reset_index(drop=True)


def _deduplicate_stage2(df: pd.DataFrame) -> pd.DataFrame:
    """Second deduplication stage on three columns."""

    subset = ["id", "target_chembl_id", "name"]
    return df.drop_duplicates(subset=subset, keep="first").reset_index(drop=True)


def _deduplicate_final(df: pd.DataFrame) -> pd.DataFrame:
    """Final deduplication stage on two columns."""

    subset = ["id", "name"]
    return df.drop_duplicates(subset=subset, keep="first").reset_index(drop=True)


@dataclass(frozen=True)
class _ProcessingStats:
    combined_rows: int
    stage1_rows: int
    stage2_rows: int
    final_rows: int


def _postprocess_frame(df: pd.DataFrame) -> tuple[pd.DataFrame, _ProcessingStats]:
    """Execute the isoform transformation pipeline."""

    projected = df.loc[:, TARGET_INPUT_COLUMNS].copy()
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

    exploded = projected.explode("triples", ignore_index=True)

    if exploded.empty:
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
        triples_df = pd.DataFrame(exploded["triples"].tolist())
        expanded = pd.concat(
            [
                exploded[["uniprot_id_primary", "target_chembl_id"]].reset_index(
                    drop=True
                ),
                triples_df,
            ],
            axis=1,
        )

    expanded = expanded.astype({col: object for col in expanded.columns})
    expanded["tokens"] = expanded["synonym"].map(_tokenize_synonym)

    names_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "name"]].copy()
    names_raw["name"] = names_raw["name"].map(lambda value: _to_text(value).strip())
    names_clean = names_raw[
        (names_raw["name"] != "")
        & (names_raw["name"] != "n/a")
        & (names_raw["name"] != "none")
    ].reset_index(drop=True)

    syns_raw = expanded[["id", "uniprot_id_primary", "target_chembl_id", "tokens"]]
    syns_exploded = syns_raw.explode("tokens", ignore_index=True)
    if syns_exploded.empty:
        syns_clean = pd.DataFrame(columns=names_clean.columns)
    else:
        syns_exploded["tokens"] = syns_exploded["tokens"].map(
            lambda value: _to_text(value).strip()
        )
        syns_filtered = syns_exploded[
            (syns_exploded["tokens"] != "")
            & (syns_exploded["tokens"] != "n/a")
            & (syns_exploded["tokens"] != "none")
        ].reset_index(drop=True)
        syns_clean = syns_filtered.rename(columns={"tokens": "name"})

    combined = pd.concat([names_clean, syns_clean], ignore_index=True)
    combined = combined.astype({col: object for col in combined.columns})
    stage1 = _deduplicate_stage1(combined)
    sorted_stage = _sort_for_stage2(stage1)
    stage2 = _deduplicate_stage2(sorted_stage)
    final = _deduplicate_final(stage2)
    final = final.loc[:, FINAL_COLUMN_ORDER]

    stats = _ProcessingStats(
        combined_rows=len(combined),
        stage1_rows=len(stage1),
        stage2_rows=len(stage2),
        final_rows=len(final),
    )
    return final, stats


def _read_input_csv(path: Path) -> pd.DataFrame:
    """Read the source CSV trying UTF encodings before CP1252."""

    last_error: Exception | None = None
    for encoding in CSV_ENCODING_PREFERENCE:
        try:
            return pd.read_csv(
                path,
                dtype=str,
                keep_default_na=False,
                na_filter=False,
                encoding=encoding,
            )
        except UnicodeDecodeError as exc:  # pragma: no cover - encoding fallback
            last_error = exc
            continue
    if last_error is not None:  # pragma: no cover - all fallbacks failed
        raise last_error
    raise FileNotFoundError(path)


def _resolve_input(path: str | Path | None) -> Path:
    if path is not None:
        return Path(path)

    directory = DEFAULT_SEARCH_DIR
    if not directory.exists():
        raise FileNotFoundError("No target output found in default directory")
    candidates = sorted(
        (candidate for candidate in directory.glob("output.target_*.csv")),
        key=lambda item: item.stat().st_mtime,
    )
    if not candidates:
        raise FileNotFoundError("No output.target_*.csv files located")
    return candidates[-1]


def process_targets(
    input_csv: str | None = None,
    output_csv: str | None = None,
    verbose: bool = True,
) -> Path:
    """Post-process the target isoform table.

    Parameters
    ----------
    input_csv : str, optional
        Source CSV produced by the target pipeline. When ``None`` the latest
        ``output.target_*.csv`` inside :data:`DEFAULT_SEARCH_DIR` is used.
    output_csv : str, optional
        Destination CSV path. Defaults to ``isoform.<input name>``.
    verbose : bool, optional
        Emit a concise progress summary to stdout when ``True``.
    """

    input_path = _resolve_input(input_csv)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    output_path = Path(output_csv) if output_csv else input_path.with_name(
        f"isoform.{input_path.name}"
    )

    if verbose:
        print(f"[isoform] reading {input_path}")

    frame = _read_input_csv(input_path)
    missing = [col for col in TARGET_INPUT_COLUMNS if col not in frame.columns]
    if missing:
        missing_fmt = ", ".join(missing)
        raise KeyError(f"Missing required columns: {missing_fmt}")

    processed, stats = _postprocess_frame(frame)

    if verbose:
        print(
            "[isoform] rows before dedup: {before}".format(
                before=stats.combined_rows
            )
        )
        print(
            "[isoform] stage1 -> {rows}".format(
                rows=stats.stage1_rows,
            )
        )
        print(
            "[isoform] stage2 -> {rows}".format(
                rows=stats.stage2_rows,
            )
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    processed.to_csv(output_path, index=False, encoding="utf-8")

    if verbose:
        print(
            "[isoform] wrote {rows} rows to {path}".format(
                rows=stats.final_rows, path=output_path
            )
        )

    return output_path

