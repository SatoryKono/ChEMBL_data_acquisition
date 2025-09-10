"""Utilities for document type classification.

Provides :func:`classify_dataframe` mirroring the behaviour expected by the
unit tests. It leverages :mod:`library.document_type_classifier` for the actual
scoring logic.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import argparse
import logging
from pathlib import Path

import pandas as pd
from library.config import Config, ensure_dirs
from library.cli import apply_config_overrides, configure_logging
from library import io

from library.document_type_classifier import compute_scores, decide_label

logger = logging.getLogger(__name__)


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


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for document type classification."""


    parser = base_parser(__doc__ or "Document type classification", column="chembl_id")
    args = parser.parse_args(argv)
    cfg: Config = apply_config_overrides(args, parser, args.config)
    ensure_dirs(cfg)
    configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)


    df_in = pd.read_csv(args.input_csv, sep=args.sep, encoding=args.encoding)
    df_out = classify_dataframe(df_in)
    output = args.output_csv or io.default_output_path(args.input_csv)
    df_out.to_csv(output, index=False, sep=args.sep, encoding=args.encoding)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
