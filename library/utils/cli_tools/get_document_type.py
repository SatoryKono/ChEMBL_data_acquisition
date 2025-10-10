"""Utilities for document type classification.

Provides :func:`classify_dataframe` mirroring the behaviour expected by the
unit tests. It leverages :mod:`library.pipelines.document.type_classifier` for the actual
scoring logic.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import pandas as pd

from library import cli, io
from library.cli import (
    build_parser as base_parser,
)
from library.cli import (
    configure_logger,
)
from library.common.log import logger
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.pipelines.document.type_classifier import compute_scores, decide_label


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
    """Command-line entry point for document type classification.

    Notes
    -----
    Relative paths honour ``--base-path``, ``--input-dir`` and ``--output-dir``.
    """

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
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of rows to process",
    )
    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be a positive integer")
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        log_cfg.run_id = run_id_value
    log_cfg.level = args.log_level
    logger_inst = configure_logger(log_cfg)
    logger_inst.info("pipeline_start", run_id=log_cfg.run_id)

    try:
        cfg: Config = cli.apply_config_overrides(
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
                "limit": "doc_type.limit",
            },
        )
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger_inst.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger_inst.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger_inst = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger_inst.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error(
            "directory_setup_failed",
            error=str(exc),
            output=str(args.output_csv),
        )
        logger_inst.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    limit_cfg = cfg.doc_type.limit
    if limit_cfg is not None and limit_cfg < 0:
        logger.error(
            "invalid_limit",
            section="doc_type.limit",
            limit=limit_cfg,
        )
        logger_inst.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    df_in = pd.read_csv(
        args.input_csv,
        sep=args.sep,
        encoding=args.encoding,
        nrows=limit_cfg,
    )
    if limit_cfg is not None:
        logger_inst.info("process_limit", limit=len(df_in))
    df_out = classify_dataframe(
        df_in,
        weights=cfg.doc_type.weights,
        thresholds=cfg.doc_type.thresholds,
    )
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df_out.to_csv(output, index=False, sep=args.sep, encoding=args.encoding)
    logger_inst.info("pipeline_done", run_id=log_cfg.run_id)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
