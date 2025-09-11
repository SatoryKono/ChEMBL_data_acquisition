"""CLI for combining ChEMBL initialisation Excel sources.

Exports pair tables without merging:

- ``pairs_same_document.csv`` from sheet ``pairs_same_doc`` in ``--same-doc``
- ``pairs_independent.csv`` and ``pairs_non_independent.csv`` derived from
  ``step5_pairs`` in ``--all-doc`` based on the ``INDEPENDENT`` flag.

Additionally, for each pair table the corresponding ``activity``, ``assay``,
``document``, ``target`` and ``testitem`` entries are exported with matching
suffixes, for example ``activity_independent.csv`` or
``assay_same_document.csv``.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library import input_initialisation_library as lib  # noqa: E402
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
from library.table_quality import analyze_table_quality  # noqa: E402


def run(cfg: Config, args: argparse.Namespace) -> int:
    """Execute table combination routine.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command line arguments.

    Returns
    -------
    int
        Zero on success, non-zero otherwise.

    """
    try:
        # ``same_doc`` and ``all_doc`` may be ``None`` when the caller bypasses
        # :func:`main` and invokes :func:`run` directly.  Guard against this to
        # avoid ``AttributeError`` when accessing ``Path.exists``.
        if args.same_doc is None or not args.same_doc.exists():
            raise FileNotFoundError(f"{args.same_doc} does not exist")
        if args.all_doc is None or not args.all_doc.exists():
            raise FileNotFoundError(f"{args.all_doc} does not exist")

        out_dir = Path(args.out_dir or cfg.init.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Loading source workbooks")
        same = lib.load_same_doc(args.same_doc)
        all_ = lib.load_all_doc(args.all_doc)

        logger.info("Combining tables")
        tables = lib.build_combined_tables(
            same,
            all_,
            dictionary_dir=args.dictionary_dir,
            status_csv=cfg.resources.status_csv,
            targets_type_csv=cfg.resources.targets_type_csv,
        )
        logger.info("Generating entity tables for pair segments")
        tables = lib.generate_pair_entity_tables(
            tables,
            {
                "pairs_independent": "independent",
                "pairs_non_independent": "non_independent",
                "pairs_same_document": "same_document",
            },
        )
        logger.info("Computing status percentages")
        for key, df in list(tables.items()):
            if not key.endswith("_status"):
                continue
            if "Filtered.new" not in df.columns:
                logger.warning("table '%s' lacks Filtered.new; skipping", key)
                continue
            entity = key.split("_")[0]
            tables[key] = lib.compute_status_statistics(df, entity)

        logger.info("Saving output")
        paths = lib.save_tables(tables, out_dir, cfg, fmt=args.format)
        # Ensure that files were actually written to disk
        missing = [str(p) for p in paths.values() if not p.exists()]
        if missing:
            raise RuntimeError("failed to write output files: " + ", ".join(missing))

        logger.info("Generating data quality reports")
        report_dir = out_dir / "data_validity_report"
        report_dir.mkdir(parents=True, exist_ok=True)
        for entity, path in paths.items():
            logger.info("Profiling %s", entity)
            analyze_table_quality(path, table_name=str(report_dir / path.stem))

        logger.info("Saved %d tables and quality reports to %s", len(paths), out_dir)
        return 0
    except KeyError as exc:
        logger.error("required table '%s' missing", exc.args[0])
        return 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create argument parser."""
    parser, log_cfg = base_parser(__doc__ or "Input initialisation", column="chembl_id")
    parser.add_argument(
        "--same-doc",
        type=Path,
        help="Path to same document workbook (default: config init.same_doc)",
    )
    parser.add_argument(
        "--all-doc",
        type=Path,
        help="Path to all document workbook (default: config init.all_doc)",
    )
    parser.add_argument(
        "--dictionary-dir",
        type=Path,
        default=None,
        help=(
            "Directory with targets_type.csv, citation_fraction.csv and status.csv "
            "(default: config resources.dictionary_dir)"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory (default: config init.output_dir)",
    )
    parser.add_argument(
        "--format", choices=["csv"], default="csv", help="Output format"
    )
    parser.set_defaults(func=run)
    return parser, log_cfg


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()

    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    try:
        cfg: Config = apply_config_overrides(
            args,
            parser,
            args.config,
            mapping={
                "same_doc": "init.same_doc",
                "all_doc": "init.all_doc",
                "out_dir": "init.output_dir",
                "dictionary_dir": "resources.dictionary_dir",
            },
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
            return 0
        ensure_dirs(cfg)

        # ``dictionary_dir`` is optional and may be undefined in the
        # configuration.  Convert provided paths to :class:`Path` objects only
        # when values are available to avoid ``TypeError`` on ``None``.
        args.same_doc = Path(args.same_doc) if args.same_doc is not None else None
        args.all_doc = Path(args.all_doc) if args.all_doc is not None else None
        # ``out_dir`` is mandatory but may originate from the configuration
        # file.  Resolve to :class:`Path` lazily to avoid ``TypeError`` when both
        # sources are missing.
        args.out_dir = (
            Path(args.out_dir)
            if args.out_dir is not None
            else Path(cfg.init.output_dir)
        )
        args.dictionary_dir = (
            Path(args.dictionary_dir) if args.dictionary_dir is not None else None
        )
        if args.same_doc is None or args.all_doc is None:
            raise ValueError("same_doc and all_doc paths must be provided")

        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    exit_code = int(args.func(cfg, args))
    if exit_code == 0:
        logger.info("pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"})
    else:
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
    return exit_code


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
