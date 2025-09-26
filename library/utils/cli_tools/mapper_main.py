"""Command line interface for mapping ChEMBL target IDs to UniProt accessions."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from urllib.error import URLError

import pandas as pd

from library import io
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    configure_logger,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import Config, ensure_dirs, print_config
from library.log import logger
from library.mapper_library import map_chembl_to_uniprot


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
            logger.error("read_fail", error=str(exc))
            return 1
        if args.column not in df.columns:
            logger.error("missing_column", column=args.column, path=str(args.input_csv))
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
                    logger.info(
                        "mapped",
                        chembl_id=str(chembl_id),
                        uniprot_id=uniprot_id,
                    )
                else:
                    logger.warning("uniprot_id_missing", chembl_id=str(chembl_id))
            except (ValueError, TimeoutError, URLError) as exc:
                logger.warning("map_failed", chembl_id=str(chembl_id), error=str(exc))
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
            logger.info("write_done", path=str(output))
        except OSError as exc:
            logger.error("write_fail", error=str(exc))
            return 1
        return 0
    except Exception as exc:  # pragma: no cover - defensive
        logger.exception("run_fail", exc=exc)
        return 1


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
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = apply_config_overrides(args, parser, args.config)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("config_error", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("setup_fail", error=str(exc))
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1
    exit_code: int = args.func(cfg, args)
    if exit_code == 0:
        logger.info("pipeline_done", run_id=log_cfg.run_id)
    else:
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
    return exit_code


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
