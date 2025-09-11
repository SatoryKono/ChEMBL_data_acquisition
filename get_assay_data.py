"""Command line interface for retrieving ChEMBL assay data."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import requests
from pandera.errors import SchemaErrors
from schemas.assays import AssaysSchema
from library.normalization import normalize_assays
from library.sidecar import SidecarErrors

from library import assay_postprocessing as ap
from library import chembl_library as cl
from library import io
from library.cli import build_parser as base_parser, configure_logging
from library.config import load_config
from library.table_quality import analyze_table_quality

logger = logging.getLogger(__name__)


def run_chembl(args: argparse.Namespace) -> int:
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
            sep=args.sep,
            encoding=args.encoding,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    try:
        df = cl.get_assays(ids, chunk_size=args.chunk_size)
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve assays: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv)
    df = ap.postprocess_assays(df)
    df = normalize_assays(df)
    sidecar = SidecarErrors()
    try:
        df = AssaysSchema.validate(df, lazy=True)
    except SchemaErrors as err:
        err.failure_cases.to_csv(output.with_name("failure_cases.csv"), index=False)
        for row in err.failure_cases.to_dict(orient="records"):
            sidecar.add_error(row)
        bad_idx = err.failure_cases["index"].dropna().unique()
        logger.warning("schema validation failed for %d rows", len(bad_idx))
        df = df.drop(index=bad_idx)
    sidecar.save(output.with_suffix(".errors.csv"))
    try:
        io.write_csv(df, output, sep=args.sep, encoding=args.encoding)
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
    args = parser.parse_args(argv)
    config = load_config(args.config)
    if args.output_csv is None:
        out_dir = config.get("output", {}).get("data_dir")
        if out_dir:
            args.output_csv = (
                Path(out_dir) / io.default_output_path(args.input_csv).name
            )
    configure_logging(args.log_level)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
