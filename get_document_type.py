"""Utilities for document type classification.

Provides :func:`classify_dataframe` mirroring the behaviour expected by the
unit tests. It leverages :mod:`library.document_type_classifier` for the actual
scoring logic.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
from library.config import Config, load_config
from library.cli import configure_logging

from library.document_type_classifier import compute_scores, decide_label


cfg: Config = load_config()


def _split_terms(value: object) -> Iterable[str]:
    """Split a semicolon separated string into terms."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    return [t.strip().lower() for t in str(value).split(";") if t]


def classify_dataframe(
    df: pd.DataFrame,
    *,
    min_review_score: int = 1,
    min_unknown_score: int = 2,
    min_experimental_score: int = 1,
) -> pd.DataFrame:
    """Classify rows in ``df`` into document types.

    Parameters
    ----------
    df:
        Input dataframe containing publication type columns.
    min_review_score, min_unknown_score, min_experimental_score:
        Thresholds for :func:`decide_label`.

    Returns
    -------
    pandas.DataFrame
        Original dataframe with an added ``class_label`` column.

    """
    result = df.copy()

    def _classify(row: pd.Series) -> str:
        scores = compute_scores(
            _split_terms(row.get("PubMed.PublicationType")),
            _split_terms(row.get("scholar.PublicationTypes")),
            _split_terms(row.get("OpenAlex.PublicationTypes")),
        )
        return decide_label(
            scores,
            min_review_score=min_review_score,
            min_unknown_score=min_unknown_score,
            min_experimental_score=min_experimental_score,
        )

    result["class_label"] = df.apply(_classify, axis=1)
    return result


def main() -> int:  # pragma: no cover - simple CLI
    """Command-line entry point for document type classification.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument("--input", type=Path, required=True, help="Input CSV")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()
    global cfg
    cfg = load_config(args.config)
    if args.log_level == "INFO":
        args.log_level = cfg.log.level
    configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)

    df_in = pd.read_csv(args.input, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding)
    df_out = classify_dataframe(df_in)
    df_out.to_csv(
        args.output, index=False, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
