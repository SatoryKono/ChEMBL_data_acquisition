"""Target post-processing helpers.

This module implements a lightweight transformation that mirrors the
Power Query workbook used during the early iterations of the pipeline.
The intent is to keep a deterministic, fully scripted variant that can
be exercised in automated tests without relying on Excel.
"""

from __future__ import annotations

import logging
import re
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd

LOGGER = logging.getLogger(__name__)

# Power Query relied on a UTF-8 CSV exported from Excel, but there are a few
# legacy snapshots encoded in various Windows code pages.  The function below
# mimics the behaviour by attempting a handful of common encodings before
# giving up.
DEFAULT_ENCODINGS: Sequence[str] = (
    "utf-8-sig",
    "utf-8",
    "cp1251",
    "cp1252",
    "latin-1",
)

# Tokens are extracted with a simple alphanumeric matcher.  Hyphens and other
# separators are treated as token boundaries which mirrors the SynExpand step
# present in the original Power Query file.
TOKEN_PATTERN = re.compile(r"[A-Za-z0-9']+")


def _normalise_series(series: pd.Series, *, lowercase: bool) -> pd.Series:
    """Return a cleaned copy of ``series``.

    Empty values and sentinel strings such as ``"None"`` or ``"n/a"`` are
    converted to the empty string.  When ``lowercase`` is true, the
    normalised values are converted to lower case to replicate the synonym
    handling that Power Query performed.
    """

    if series.empty:
        return series.copy()

    cleaned = (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.replace(r"\s+", " ", regex=True)
    )
    cleaned = cleaned.replace(
        to_replace=r"^(?i)(?:na|n/a|none|n\\.?a\\.?|n\s*/\s*a)$",
        value="",
        regex=True,
    )
    if lowercase:
        cleaned = cleaned.str.lower()
    return cleaned


def _split_pipe(value: str) -> list[str]:
    if not value:
        return []
    return [part.strip() for part in value.split("|") if part.strip()]


def _split_synonyms(value: str) -> list[str]:
    if not value:
        return []
    parts = [segment.strip() for segment in value.split(":") if segment.strip()]
    unique: list[str] = []
    for part in parts:
        lower = part.lower()
        if lower in {"n/a", "none", "na"}:
            continue
        if part not in unique:
            unique.append(part)
    return unique


def _syn_expand(term: str) -> list[str]:
    if not term:
        return []
    tokens = TOKEN_PATTERN.findall(term.lower())
    # Preserve order whilst removing duplicates.
    return list(dict.fromkeys(tokens))


def _read_csv(path: Path, *, encodings: Iterable[str], sep: str) -> pd.DataFrame:
    last_error: Exception | None = None
    for encoding in encodings:
        try:
            LOGGER.info("Reading target isoform table from %s (encoding=%s)", path, encoding)
            return pd.read_csv(path, dtype=str, encoding=encoding, sep=sep)
        except UnicodeDecodeError as exc:  # pragma: no cover - defensive
            last_error = exc
            continue
    if last_error is not None:
        raise last_error
    raise UnicodeDecodeError("utf-8", b"", 0, 1, "unable to decode input")


def process_targets(
    source: Path | str,
    *,
    output_dir: Path | str | None = None,
    sep: str = ",",
    encodings: Sequence[str] | None = None,
) -> Path:
    """Process the aggregated target CSV produced by :mod:`scripts.get_target_data`.

    The transformation mirrors the legacy Power Query workbook.  It reads the
    CSV, normalises name and synonym columns, expands pipe-separated fields
    into per-isoform entries, performs a synonym tokenisation step and writes
    the result to an ``isoform.output.<original>`` file located alongside the
    source file (or in ``output_dir`` when provided).

    Parameters
    ----------
    source:
        Path to the CSV file that should be processed.
    output_dir:
        Optional destination directory.  When omitted the output is written to
        the same directory as ``source``.
    sep:
        Field separator used in the CSV.  Defaults to a comma.
    encodings:
        Ordered list of encodings that will be attempted while reading the
        source file.  The first successful attempt stops the search.

    Returns
    -------
    pathlib.Path
        Path to the generated output file.
    """

    input_path = Path(source)
    if not input_path.exists():  # pragma: no cover - guarded by callers
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    if encodings is None:
        encodings = DEFAULT_ENCODINGS

    df = _read_csv(input_path, encodings=encodings, sep=sep)

    name_columns = [col for col in df.columns if col.endswith("_names")]
    synonym_columns = [col for col in df.columns if col.endswith("_synonyms")]

    for col in name_columns:
        df[col] = _normalise_series(df[col], lowercase=False)
    for col in synonym_columns:
        df[col] = _normalise_series(df[col], lowercase=True)

    records: list[dict[str, str]] = []
    for row in df.itertuples(index=False):
        row_dict = row._asdict()
        target_id = row_dict.get("target_chembl_id", "")
        isoform_ids = _split_pipe(row_dict.get("isoform_ids", ""))
        isoform_names = _split_pipe(row_dict.get("isoform_names", ""))
        isoform_synonyms = [
            _split_synonyms(value)
            for value in _split_pipe(row_dict.get("isoform_synonyms", ""))
        ]
        for iso_id, iso_name, iso_syns in zip_longest(
            isoform_ids, isoform_names, isoform_synonyms, fillvalue=""
        ):
            # Normalise placeholders that occasionally sneak into the CSV.
            iso_id = iso_id.strip() if isinstance(iso_id, str) else ""
            if iso_id.lower() in {"n/a", "none", "na"}:
                iso_id = ""
            iso_name = iso_name.strip() if isinstance(iso_name, str) else ""
            if iso_name.lower() in {"n/a", "none", "na"}:
                iso_name = ""

            if not iso_name and not iso_syns:
                continue

            candidates = list(dict.fromkeys([iso_name] if iso_name else []))
            for synonym in iso_syns or []:
                if synonym and synonym not in candidates:
                    candidates.append(synonym)

            for term in candidates:
                tokens = _syn_expand(term)
                if not tokens:
                    continue
                for token in tokens:
                    if not token:
                        continue
                    records.append(
                        {
                            "target_chembl_id": target_id,
                            "isoform_id": iso_id,
                            "isoform_name": iso_name,
                            "term": term,
                            "token": token,
                        }
                    )

    if records:
        result = pd.DataFrame.from_records(records)
    else:
        result = pd.DataFrame(
            columns=["target_chembl_id", "isoform_id", "isoform_name", "term", "token"]
        )

    initial_size = len(result)
    if initial_size:
        result = result.drop_duplicates(
            subset=["target_chembl_id", "isoform_id", "term", "token"],
            keep="first",
            ignore_index=True,
        )
    LOGGER.info("Isoform synonym tokens: %d -> %d rows", initial_size, len(result))

    if len(result):
        result = result.sort_values(
            by=["target_chembl_id", "isoform_id", "term", "token"],
            kind="mergesort",
            ignore_index=True,
        )

    output_directory = Path(output_dir) if output_dir is not None else input_path.parent
    output_path = output_directory / f"isoform.output.{input_path.name}"
    LOGGER.info("Writing processed target isoforms to %s", output_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False, encoding="utf-8", line_terminator="\n")
    return output_path
