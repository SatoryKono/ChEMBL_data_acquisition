"""Command line interface for mapping ChEMBL target IDs to UniProt accessions."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from urllib.error import URLError

import pandas as pd

from library import cli
from library import io
from library.cli import (
    LoggerConfig,
    configure_logger,
    positive_int,
)
from library.cli import (
    build_parser as base_parser,
)
from library.config import Config, ensure_dirs, print_config
from library.log import logger
from library.mapper_batch_library import map_chembl_ids_to_uniprot


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

        def _normalize(value: object) -> str | None:
            if pd.isna(value):
                return None
            text = str(value).strip()
            return text or None

        ids_to_map: list[str] = []
        for value in df[args.column]:
            normalized = _normalize(value)
            if normalized is not None:
                ids_to_map.append(normalized)

        try:
            mappings = map_chembl_ids_to_uniprot(
                ids_to_map,
                cfg.uniprot_mapping,
                batch_size=args.chunk_size,
                rps=args.rps,
                max_workers=args.workers,
            )
        except (ValueError, TimeoutError, URLError) as exc:
            logger.warning("map_failed", error=str(exc))
            df["mapping_uniprot_id"] = [None for _ in df[args.column]]
        else:
            for chembl_id in ids_to_map:
                uniprot_id = mappings.get(chembl_id)
                if uniprot_id:
                    logger.info(
                        "mapped",
                        chembl_id=chembl_id,
                        uniprot_id=uniprot_id,
                    )
                else:
                    logger.warning("uniprot_id_missing", chembl_id=chembl_id)
            mapped_values: list[str | None] = []
            for value in df[args.column]:
                normalized = _normalize(value)
                mapped_values.append(
                    mappings.get(normalized) if normalized is not None else None
                )
            df["mapping_uniprot_id"] = mapped_values
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
    parser.add_argument("--rps", type=float, default=1.0, help="Max requests per second")
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="Number of worker threads",
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
        cfg: Config = cli.apply_config_overrides(args, parser, args.config)
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg)
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
