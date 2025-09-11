"""Utilities for document type classification.

Provides :func:`classify_dataframe` mirroring the behaviour expected by the
unit tests. It leverages :mod:`library.document_type_classifier` for the actual
scoring logic.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import logging

import pandas as pd
from library.config import Config, ensure_dirs, print_config
from library.cli import (
    apply_config_overrides,
    build_parser as base_parser,
    configure_logger,
)
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
    weights: Mapping[str, int] | None = None,
    thresholds: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Classify rows in ``df`` into document types.

    Parameters
    ----------
    df:
        Input dataframe containing publication type columns.
    weights:
        Optional mapping of source names to weights for
        :func:`compute_scores`.
    thresholds:
        Optional mapping with keys ``review``, ``experimental`` and
        ``unknown`` specifying minimum scores for :func:`decide_label`.

    Returns
    -------
    pandas.DataFrame
        Original dataframe with an added ``class_label`` column.

    """
    result = df.copy()
    thresh = thresholds or {"review": 1, "experimental": 1, "unknown": 2}

    def _classify(row: pd.Series) -> str:
        scores = compute_scores(
            _split_terms(row.get("PubMed.PublicationType")),
            _split_terms(row.get("scholar.PublicationTypes")),
            _split_terms(row.get("OpenAlex.PublicationTypes")),
            weights=weights,
        )
        return decide_label(
            scores,
            min_review_score=thresh.get("review", 1),
            min_unknown_score=thresh.get("unknown", 2),
            min_experimental_score=thresh.get("experimental", 1),
        )

    result["class_label"] = df.apply(_classify, axis=1)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point for document type classification."""

    parser, log_cfg = base_parser(
        __doc__ or "Document type classification", column="chembl_id"
    )
    parser.add_argument(
        "--weight-pubmed",
        dest="weight_pubmed",
        type=int,
        default=None,
        help="Weight for PubMed source",
    )
    parser.add_argument(
        "--weight-scholar",
        dest="weight_scholar",
        type=int,
        default=None,
        help="Weight for Scholar source",
    )
    parser.add_argument(
        "--weight-openalex",
        dest="weight_openalex",
        type=int,
        default=None,
        help="Weight for OpenAlex source",
    )
    parser.add_argument(
        "--threshold-review",
        dest="threshold_review",
        type=int,
        default=None,
        help="Minimum score for review label",
    )
    parser.add_argument(
        "--threshold-experimental",
        dest="threshold_experimental",
        type=int,
        default=None,
        help="Minimum score for experimental label",
    )
    parser.add_argument(
        "--threshold-unknown",
        dest="threshold_unknown",
        type=int,
        default=None,
        help="Minimum score for unknown label",
    )
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger_inst = configure_logger(log_cfg)
    logger_inst.info(
        "pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"}
    )

    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "weight_pubmed": "doc_type.weights.pubmed",
                "weight_scholar": "doc_type.weights.scholar",
                "weight_openalex": "doc_type.weights.openalex",
                "threshold_review": "doc_type.thresholds.review",
                "threshold_experimental": "doc_type.thresholds.experimental",
                "threshold_unknown": "doc_type.thresholds.unknown",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger_inst.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
            return 0
        ensure_dirs(cfg)
        logger_inst = configure_logger(
            log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt
        )
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger_inst.info(
            "pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"}
        )
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger_inst.info(
            "pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"}
        )
        return 1

    df_in = pd.read_csv(args.input_csv, sep=args.sep, encoding=args.encoding)
    df_out = classify_dataframe(
        df_in,
        weights=cfg.doc_type.weights,
        thresholds=cfg.doc_type.thresholds,
    )
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df_out.to_csv(output, index=False, sep=args.sep, encoding=args.encoding)
    logger_inst.info("pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"})
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
