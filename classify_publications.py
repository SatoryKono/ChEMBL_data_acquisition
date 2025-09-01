"""CLI tool for classifying scientific publications.

The script reads a CSV file with publication metadata, assigns each record to
``review``, ``experimental`` or ``unknown`` class and writes results to a new
CSV file with additional debug columns.
"""
from __future__ import annotations

import argparse
import logging
import random
from typing import Sequence

import pandas as pd

from doc_classifier.classifier import compute_scores, decide_label
from doc_classifier.terms import parse_terms

REQUIRED_COLUMNS = [
    "title",
    "abstract",
    "PubMed.PublicationType",
    "scholar.PublicationTypes",
    "OpenAlex.PublicationTypes",
    "OpenAlex.TypeCrossref",
]


def _is_empty(series: pd.Series) -> pd.Series:
    """Return boolean Series indicating empty strings or NaNs."""
    return series.isna() | (series.astype(str).str.strip() == "")


def _read_csv_auto(path: str, encodings: Sequence[str], seps: Sequence[str]) -> tuple[pd.DataFrame, str, str]:
    """Read a CSV file trying multiple encodings and separators."""
    last_error: Exception | None = None
    for enc in encodings:
        for sep in seps:
            try:
                df = pd.read_csv(path, sep=sep, encoding=enc)
            except Exception as err:  # pragma: no cover - diagnostic path
                last_error = err
                continue
            if set(REQUIRED_COLUMNS).issubset(df.columns):
                return df, enc, sep
    raise ValueError(f"Cannot read input file with provided encodings/separators: {last_error}")


def classify_dataframe(
    df: pd.DataFrame,
    *,
    min_review_score: int,
    min_unknown_score: int,
    min_experimental_score: int,
) -> pd.DataFrame:
    """Classify DataFrame rows and append debug columns."""
    df = df.copy()

    pubmed_terms = df["PubMed.PublicationType"].map(parse_terms)
    scholar_terms = df["scholar.PublicationTypes"].map(parse_terms)
    oa_pt = df["OpenAlex.PublicationTypes"].map(parse_terms)
    oa_cr = df["OpenAlex.TypeCrossref"].map(parse_terms)
    openalex_terms = oa_pt.combine(oa_cr, lambda a, b: sorted(set(a) | set(b)))

    df["debug.pubmed_terms"] = pubmed_terms.map(lambda t: "|".join(t))
    df["debug.scholar_terms"] = scholar_terms.map(lambda t: "|".join(t))
    df["debug.openalex_terms"] = openalex_terms.map(lambda t: "|".join(t))

    def classify_row(row):
        scores = compute_scores(
            pubmed_terms[row.name], scholar_terms[row.name], openalex_terms[row.name]
        )
        label = decide_label(
            scores,
            min_review_score=min_review_score,
            min_unknown_score=min_unknown_score,
            min_experimental_score=min_experimental_score,
        )
        return pd.Series(
            {
                "class_label": label,
                "debug.scores.review": scores["review"],
                "debug.scores.experimental": scores["experimental"],
                "debug.scores.unknown": scores["unknown"],
            }
        )

    results = df.apply(classify_row, axis=1)
    df = pd.concat([df, results], axis=1)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="Classify publication records")
    parser.add_argument("--input", required=True, help="Path to input CSV file")
    parser.add_argument("--output", required=True, help="Path to output CSV file")
    parser.add_argument("--min-review-score", type=int, default=1)
    parser.add_argument("--min-unknown-score", type=int, default=2)
    parser.add_argument("--min-experimental-score", type=int, default=1)
    parser.add_argument("--random-seed", type=int, default=0)
    parser.add_argument("--encoding-fallbacks", default="utf-8-sig,cp1251,latin1")
    parser.add_argument("--delimiters", default=",;|\t")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    random.seed(args.random_seed)

    encodings = [e.strip() for e in args.encoding_fallbacks.split(",") if e.strip()]
    seps = [s for s in args.delimiters]
    df, encoding, sep = _read_csv_auto(args.input, encodings, seps)

    logging.info("Detected encoding %s and delimiter '%s'", encoding, sep)

    for col in REQUIRED_COLUMNS:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' not found in input")

    classified = classify_dataframe(
        df,
        min_review_score=args.min_review_score,
        min_unknown_score=args.min_unknown_score,
        min_experimental_score=args.min_experimental_score,
    )

    counts = classified["class_label"].value_counts().sort_index()
    total = len(classified)
    for label, cnt in counts.items():
        print(f"{label}: {cnt} ({cnt / total:.1%})")

    pubmed_empty = _is_empty(df["PubMed.PublicationType"]).mean() * 100
    scholar_empty = _is_empty(df["scholar.PublicationTypes"]).mean() * 100
    openalex_empty = (
        _is_empty(df["OpenAlex.PublicationTypes"]) & _is_empty(df["OpenAlex.TypeCrossref"])
    ).mean() * 100
    zero_signals = (
        (classified[["debug.scores.review", "debug.scores.experimental", "debug.scores.unknown"]].sum(axis=1) == 0)
    ).sum()

    logging.info(
        "Empty fields - PubMed: %.1f%%, Scholar: %.1f%%, OpenAlex: %.1f%%",
        pubmed_empty,
        scholar_empty,
        openalex_empty,
    )
    logging.info("Rows with zero signals: %d", zero_signals)

    classified.to_csv(args.output, index=False, sep=sep, encoding=encoding)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    main()
