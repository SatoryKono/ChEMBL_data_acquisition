"""Command line interface for mapping ChEMBL target IDs to UniProt accessions."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from collections.abc import Sequence
from typing import Any
from pathlib import Path
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
from library.config import Config, ConfigError, ensure_dirs, print_config
from library.common.log import logger
from library.integration.mapper_batch_library import map_chembl_ids_to_uniprot


SUMMARY_SAMPLE_LIMIT = 5


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
        marker_values = list(cfg.io.na_markers or ())
        keep_markers = bool(getattr(cfg.io, "keep_na_markers", False))
        na_values = marker_values if marker_values and not keep_markers else None
        try:
            df = io.read_csv(
                args.input_csv,
                cfg=cfg.io,
                sep=args.sep,
                encoding=args.encoding,
                dtype=str,
                na_values=na_values,
            )
        except (FileNotFoundError, OSError) as exc:
            logger.error("read_fail", error=str(exc))
            return 1
        if args.column not in df.columns:
            logger.error("missing_column", column=args.column, path=str(args.input_csv))
            return 1

        marker_set = set(marker_values)

        def _normalize(value: object) -> str | None:
            if pd.isna(value):
                return None
            text = str(value).strip()
            if not text:
                return None
            if not keep_markers and text in marker_set:
                return None
            return text

        row_ids: list[str | None] = []
        unique_ids: OrderedDict[str, None] = OrderedDict()
        for value in df[args.column]:
            normalized = _normalize(value)
            row_ids.append(normalized)
            if normalized is not None:
                unique_ids.setdefault(normalized, None)

        ids_to_map = list(unique_ids.keys())
        normalized_ids = [chembl_id for chembl_id in row_ids if chembl_id is not None]

        log_each = bool(getattr(args, "log_each", False))
        total_ids = len(normalized_ids)
        mapped_count = 0
        missing_count = 0
        missing_sample: list[str] = []
        mapping_failed = False

        mapping_failed = False
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
            missing_count = total_ids
            missing_sample = normalized_ids[:SUMMARY_SAMPLE_LIMIT]

            mapping_failed = True
        else:
            for chembl_id in row_ids:
                if chembl_id is None:
                    continue
                uniprot_id = mappings.get(chembl_id)
                if uniprot_id:
                    mapped_count += 1
                    if log_each:
                        logger.info(
                            "mapped",
                            chembl_id=chembl_id,
                            uniprot_id=uniprot_id,
                        )
                    else:
                        logger.debug(
                            "mapped",
                            chembl_id=chembl_id,
                            uniprot_id=uniprot_id,
                        )
                else:
                    missing_count += 1
                    if len(missing_sample) < SUMMARY_SAMPLE_LIMIT:
                        missing_sample.append(chembl_id)
                    if log_each:
                        logger.warning("uniprot_id_missing", chembl_id=chembl_id)
                    else:
                        logger.debug("uniprot_id_missing", chembl_id=chembl_id)
            mapped_values: list[str | None] = [
                mappings.get(chembl_id) if chembl_id is not None else None
                for chembl_id in row_ids
            ]
            df["mapping_uniprot_id"] = mapped_values

        summary_payload: dict[str, Any] = {
            "total": total_ids,
            "mapped": mapped_count,
            "missing": missing_count,
        }
        if missing_sample:
            summary_payload["sample_missing"] = missing_sample
        logger.info("mapper_summary", **summary_payload)
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
        return 1 if mapping_failed else 0
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
    parser.epilog = (
        "Retry failed mappings by re-running the command. "
        "Consider reducing --chunk-size or --rps when the remote service throttles requests."
    )
    parser.add_argument(
        "--key-cols",
        nargs="*",
        help="Columns used to sort the output CSV deterministically",
    )
    parser.add_argument(
        "--rps", type=float, default=1.0, help="Max requests per second"
    )
    parser.add_argument(
        "--workers",
        type=positive_int,
        default=1,
        help="Number of worker threads",
    )
    parser.add_argument(
        "--log-each",
        action="store_true",
        help="Emit per-ID mapping logs at INFO/WARN level",
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults.

    Notes
    -----
    Relative paths honour ``--base-path``, ``--input-dir`` and ``--output-dir``.
    """
    parser, log_cfg = build_parser()

    args = parser.parse_args(argv)
    input_path = getattr(args, "input_csv", None)
    output_stem = Path(input_path).stem if input_path else None
    cli.prepare_io_paths(args, output_stem=output_stem)
    run_id_value = getattr(args, "run_id", None)
    if isinstance(run_id_value, str):
        run_id_value = run_id_value.strip() or None
    if run_id_value is not None:
        log_cfg.run_id = run_id_value
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline_start", run_id=log_cfg.run_id)
    try:
        cfg: Config = cli.apply_config_overrides(args, parser, args.config)
    except (ConfigError, FileNotFoundError, ValueError) as exc:
        logger.error(
            "config_error",
            error=str(exc),
            config=str(args.config),
        )
        logger.info("pipeline_fail", run_id=log_cfg.run_id)
        return 1

    try:
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg)
            logger.info("pipeline_done", run_id=log_cfg.run_id)
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg)
    except (ValueError, TypeError) as exc:
        logger.error("config_error", error=str(exc), config=str(args.config))
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
