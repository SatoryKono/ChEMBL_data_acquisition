"""Command line interface for mapping ChEMBL target IDs to UniProt accessions."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from urllib.error import URLError

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library import io  # noqa: E402
from library.cli import (  # noqa: E402
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (  # noqa: E402
    build_parser as base_parser,
)
from library.config import Config, ensure_dirs, print_config  # noqa: E402
from library.log import logger  # noqa: E402
from library.mapper_library import map_chembl_to_uniprot  # noqa: E402


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Map ChEMBL target identifiers to UniProt accessions.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        df = io.read_csv(
            args.input_csv,
            cfg=cfg.io,
            sep=args.sep,
            encoding=args.encoding,
            dtype=str,
            na_values=["#N/A", ""],
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
            uniprot_id = map_chembl_to_uniprot(str(chembl_id), cfg.uniprot_mapping)
            uniprot_ids.append(uniprot_id)
            if uniprot_id:
                logger.info("mapped %s -> %s", chembl_id, uniprot_id)
            else:
                logger.warning("no UniProt ID for %s", chembl_id)
        except (ValueError, TimeoutError, URLError) as exc:
            logger.warning("failed to map %s: %s", chembl_id, exc)
            uniprot_ids.append(None)
    df["mapping_uniprot_id"] = uniprot_ids
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    try:
        io.write_csv(
            df,
            output,
            cfg=cfg,
            sep=args.sep,
            encoding=args.encoding,
            key_cols=args.key_cols,
        )
        logger.info("wrote %s", output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    return 0


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create the command-line argument parser."""
    parser, log_cfg = base_parser(
        "Map ChEMBL target IDs to UniProt accessions",
        column="chembl_id",
        chunk_size=1,
    )
    parser.add_argument(
        "--key-cols",
        nargs="*",
        help="Columns used to sort the output CSV deterministically",
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()

    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    try:
        cfg: Config = apply_config_overrides(args, parser, args.config)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    exit_code = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"})
    else:
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
