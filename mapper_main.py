"""Command line interface for mapping ChEMBL target IDs to UniProt accessions."""

from __future__ import annotations

import argparse
import logging
from typing import Sequence
from urllib.error import URLError

import pandas as pd
from library.config import Config, ensure_dirs

from library import io
from library.cli import (
    apply_config_overrides,
    build_parser as base_parser,
    configure_logging,
)
from library.mapper_library import map_chembl_to_uniprot

logger = logging.getLogger(__name__)


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Map ChEMBL target identifiers to UniProt accessions.

    Parameters
    ----------
    args:
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        df = io.read_csv(
            args.input_csv, cfg=cfg.io, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    if args.column not in df.columns:
        logger.error("column '%s' not found in %s", args.column, args.input_csv)
        return 1

    uniprot_ids: list[str | None] = []
    for chembl_id in df[args.column]:
        if pd.isna(chembl_id) or not str(chembl_id).strip():
            uniprot_ids.append(None)
            continue
        try:
            uniprot_id = map_chembl_to_uniprot(str(chembl_id))
            uniprot_ids.append(uniprot_id)
            if uniprot_id:
                logger.info("mapped %s -> %s", chembl_id, uniprot_id)
            else:
                logger.warning("no UniProt ID for %s", chembl_id)
        except (ValueError, TimeoutError, URLError) as exc:
            logger.warning("failed to map %s: %s", chembl_id, exc)
            uniprot_ids.append(None)
    df["mapping_uniprot_id"] = uniprot_ids
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
        logger.info("wrote %s", output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = base_parser(
        "Map ChEMBL target IDs to UniProt accessions",
        column="chembl_id",
        chunk_size=1,
    )
    parser.set_defaults(func=run)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""
    parser = build_parser()
    parser.add_argument(
        "--config", default="config.yaml", help="Path to YAML configuration file"
    )
    args = parser.parse_args(argv)
    cfg = apply_config_overrides(args, parser, args.config)
    ensure_dirs(cfg)
    configure_logging(args.log_level, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    return args.func(cfg, args)


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
