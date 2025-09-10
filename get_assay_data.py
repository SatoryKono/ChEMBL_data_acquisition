"""Command line interface for retrieving ChEMBL assay data."""

from __future__ import annotations

import argparse
import logging
from typing import Sequence

import requests
from library.config import Config, ensure_dirs, load_config

from library import assay_postprocessing as ap
from library import chembl_library as cl
from library import io
from library.cli import build_parser as base_parser, configure_logging
from library.table_quality import analyze_table_quality

logger = logging.getLogger(__name__)


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute assay retrieval from the ChEMBL API.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        ids = io.read_ids(
            args.input_csv,
            column=args.column,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        df = cl.get_assays(ids, cfg=cfg.api, chunk_size=args.chunk_size)
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve assays: %s", exc)
        return 1
    df = ap.postprocess_assays(df)
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = base_parser(
        "ChEMBL assay data utilities", column="assay_chembl_id", chunk_size=10
    )
    parser.set_defaults(func=run_chembl)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    ensure_dirs(cfg)

    default_chunk = parser.get_default("chunk_size")
    if args.chunk_size == default_chunk:
        args.chunk_size = cfg.jobs.chunk_size
    default_sep = parser.get_default("sep")
    if args.sep == default_sep:
        args.sep = cfg.io.csv_sep
    default_encoding = parser.get_default("encoding")
    if args.encoding == default_encoding:
        args.encoding = cfg.io.csv_encoding
    default_log = parser.get_default("log_level")
    if args.log_level == default_log:
        args.log_level = cfg.log.level
    configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    return args.func(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
