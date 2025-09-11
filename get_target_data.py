"""Command line interface for retrieving target data from external sources.

Example
-------
Fetch ChEMBL target information for identifiers in ``targets.csv``::

    python get_target_data.py chembl --config config.yaml --input targets.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Sequence

import pandas as pd
import requests
from library.config import (
    Config,
    ensure_dirs,
    print_config,
    _serialize_paths,
)
from library.chembl_client import ChemblClient

from library import chembl_library as cl
from library import io
from library import iuphar_library as ii
from library import target_postprocessing as tp
from library import uniprot_library as uu
from library.metadata import Stats, file_sha256, write_meta_yaml

from library.cli import (
    apply_config_overrides,
    build_root_parser,
    configure_logger,
    LoggerConfig,
)
from library.log import logger
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from pandera.errors import SchemaErrors
from schemas import TargetsSchema, normalize_targets

from library import write_csv_deterministic


def _pipe_merge(values: Sequence[str | None]) -> str:
    """Return a ``"|"``-joined string of unique, non-empty tokens.

    Parameters
    ----------
    values:
        Sequence of pipe-delimited strings to merge.

    Returns
    -------
    str
        Sorted, unique tokens separated by ``"|"``. Empty inputs yield
        an empty string.

    """
    tokens: set[str] = set()
    for value in values:
        if isinstance(value, str) and value:
            parts = [p.strip() for p in value.split("|") if p.strip()]
            tokens.update(parts)
    return "|".join(sorted(tokens))


def _first_token(value: str | None) -> str:
    """Return the first token from a pipe-delimited ``value``."""
    if isinstance(value, str) and value:
        return value.split("|")[0]
    return ""


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create and return the top-level CLI argument parser.

    The command line interface is organised into sub-commands for retrieving
    data from individual sources (UniProt, ChEMBL and IUPHAR) as well as a
    convenience ``all`` command that runs all pipelines and merges their
    outputs.
    """

    root, log_cfg = build_root_parser()
    parser = argparse.ArgumentParser(
        description="Target data utilities", parents=[root]
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ----------------------------
    # UniProt sub-command
    # ----------------------------
    uniprot = subparsers.add_parser(
        "uniprot",
        parents=[root],
        help="Extract information for UniProt accessions",
    )
    uniprot.add_argument(
        "--column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column in the input CSV containing UniProt accessions",
    )
    uniprot.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config resources.uniprot_data_dir)"
        ),
    )
    uniprot.set_defaults(func=run_uniprot)

    # ----------------------------
    # ChEMBL sub-command
    # ----------------------------
    chembl = subparsers.add_parser(
        "chembl",
        parents=[root],
        help="Retrieve target information from ChEMBL",
    )
    chembl.add_argument(
        "--column",
        default="chembl_id",
        help="Column name in the input CSV containing identifiers",
    )
    chembl.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    chembl.set_defaults(func=run_chembl)

    # ----------------------------
    # IUPHAR sub-command
    # ----------------------------
    iuphar = subparsers.add_parser(
        "iuphar",
        parents=[root],
        help="Map UniProt accessions to IUPHAR classifications",
    )
    iuphar.add_argument(
        "--target-csv",
        type=Path,
        default=None,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config resources.iuphar_target_csv)"
        ),
    )
    iuphar.add_argument(
        "--family-csv",
        type=Path,
        default=None,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config resources.iuphar_family_csv)"
        ),
    )
    iuphar.set_defaults(func=run_iuphar)

    # ----------------------------
    # Combined pipeline
    # ----------------------------
    all_cmd = subparsers.add_parser(
        "all",
        parents=[root],
        help="Run ChEMBL, UniProt and IUPHAR pipelines and merge results",
    )
    all_cmd.add_argument(
        "--chembl-out",
        dest="chembl_out",
        type=Path,
        help="Optional path to save intermediate ChEMBL data",
    )
    all_cmd.add_argument(
        "--uniprot-out",
        dest="uniprot_out",
        type=Path,
        help="Optional path to save intermediate UniProt data",
    )
    all_cmd.add_argument(
        "--iuphar-out",
        dest="iuphar_out",
        type=Path,
        help="Optional path to save intermediate IUPHAR data",
    )
    all_cmd.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing '<uniprot_id>.json' files "
            "(default: config resources.uniprot_data_dir)"
        ),
    )
    all_cmd.add_argument(
        "--target-csv",
        type=Path,
        default=None,
        help=(
            "Path to the _IUPHAR_target.csv file "
            "(default: config resources.iuphar_target_csv)"
        ),
    )
    all_cmd.add_argument(
        "--family-csv",
        type=Path,
        default=None,
        help=(
            "Path to the _IUPHAR_family.csv file "
            "(default: config resources.iuphar_family_csv)"
        ),
    )
    all_cmd.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="Timeout in seconds for each HTTP request",
    )
    all_cmd.add_argument(
        "--organism-csv",
        type=Path,
        default=None,
        help=(
            "CSV mapping 'genus' to organism 'type' for finalisation "
            "(default: config resources.organism_csv)"
        ),
    )
    all_cmd.add_argument(
        "--uniprot-column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column from ChEMBL output to use for UniProt processing",
    )
    all_cmd.set_defaults(func=run_all)

    setattr(
        parser,
        "subparsers_map",
        {
            "uniprot": uniprot,
            "chembl": chembl,
            "iuphar": iuphar,
            "all": all_cmd,
        },
    )

    return parser, log_cfg


def run_uniprot(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``uniprot`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command-line arguments specific to the ``uniprot`` sub-command.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    Tests
    -----
    The post-processing step is covered by
    :mod:`tests.test_target_postprocessing`.

    """
    try:
        df = pd.read_csv(
            args.input_csv, sep=args.sep, encoding=args.encoding, dtype=str
        )
        if args.column not in df.columns:
            raise ValueError(f"column '{args.column}' not found in {args.input_csv}")
        df = df.fillna("")
        df = df[
            (df[args.column].str.strip() != "") & (df[args.column] != "#N/A")
        ].reset_index(drop=True)
        ids = df[args.column].tolist()

        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=args.encoding, newline=""
        ) as tmp:
            writer = csv.DictWriter(tmp, fieldnames=["uniprot_id"], delimiter=args.sep)
            writer.writeheader()
            for uid in ids:
                writer.writerow({"uniprot_id": uid})
            tmp_path = Path(tmp.name)

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        data_dir = args.data_dir or cfg.resources.uniprot_data_dir
        try:
            uu.process(
                input_csv=str(tmp_path),
                output_csv=str(output),
                data_dir=data_dir,
                cfg=cfg.uniprot,
                sep=args.sep,
                encoding=args.encoding,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        out_df = pd.read_csv(output, sep=args.sep, encoding=args.encoding, dtype=str)
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(1, "mapping_uniprot_id", df["mapping_uniprot_id"].tolist())
        io.write_csv(out_df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(out_df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_chembl(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    # Set up HTTP client with proper headers and retry behaviour
    client = ChemblClient(cfg.api, cfg.retry, cfg.chembl)

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
        df = cl.get_targets(
            ids,
            cfg=cfg.api,
            mapping_cfg=cfg.uniprot_mapping,
            client=client,
            timeout=args.timeout,
        )
    except (requests.RequestException, ValueError) as exc:
        logger.error("failed to retrieve targets: %s", exc)
        return 1
    output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
    df = normalize_targets(df)
    rows_total = len(df)
    exit_code = 0
    required_cols = set(TargetsSchema.columns.keys())
    if required_cols.issubset(df.columns):
        try:
            df = TargetsSchema.validate(df, lazy=True)
        except SchemaErrors as exc:
            failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation failed; wrote %d failure cases to %s",
                len(exc.failure_cases),
                failure_path,
            )
            df = getattr(exc, "validated_data", df)
            exit_code = 1
    else:
        missing = required_cols.difference(df.columns)
        logger.warning("Skipping validation due to missing columns: %s", missing)
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = [c for c in ["target_chembl_id"] if c in df.columns]
        csv_path = write_csv_deterministic(
            df,
            output,
            col_order=[
                "target_chembl_id",
                "organism",
                "target_uniprot_id",
                "pH_dependence",
                "isoforms",
            ],
            key_cols=key_cols or None,
        )
        logger.info("Wrote %d rows to %s", rows_kept, csv_path)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1

    stats: Stats = {
        "rows_total": rows_total,
        "rows_kept": rows_kept,
        "rows_dropped": rows_dropped,
        "output_sha256": file_sha256(csv_path),
    }
    write_meta_yaml(
        csv_path=csv_path,
        command=" ".join(sys.argv),
        config_subset=_serialize_paths(cfg.to_dict()),
        inputs={"input_csv": str(args.input_csv)},
        stats=stats,
        schema="TargetsSchema",
    )
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def run_iuphar(cfg: Config, args: argparse.Namespace) -> int:
    """Execute the ``iuphar`` sub-command.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    """
    try:
        data = ii.IUPHARData.from_files(
            target_path=args.target_csv,
            family_path=args.family_csv,
            encoding=args.encoding,
        )
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        data.map_uniprot_file(
            input_path=args.input_csv,
            output_path=output,
            encoding=args.encoding,
            sep=args.sep,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(output, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_all(cfg: Config, args: argparse.Namespace) -> int:
    """Run ChEMBL, UniProt and IUPHAR pipelines and merge their outputs.

    The merged table is cleaned and normalised using
    :func:`library.target_postprocessing.postprocess_targets` followed by
    :func:`library.target_postprocessing.finalise_targets` before being
    written to disk.

    Parameters
    ----------
    cfg : Config
        Application configuration.
    args:
        Parsed command-line arguments specific to the ``all`` sub-command.

    Returns
    -------
    int
        Zero on success, non-zero on failure.

    Tests
    -----
    The post-processing step is covered by
    :mod:`tests.test_target_postprocessing`.

    """
    try:
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        chembl_out = args.chembl_out or output.with_name(output.stem + "_chembl.csv")
        uniprot_out = args.uniprot_out or output.with_name(output.stem + "_uniprot.csv")
        iuphar_out = args.iuphar_out or output.with_name(output.stem + "_iuphar.csv")

        # Run ChEMBL retrieval and capture results
        chembl_args = argparse.Namespace(
            input_csv=args.input_csv,
            output_csv=chembl_out,
            column="chembl_id",
            sep=args.sep,
            encoding=args.encoding,
            timeout=args.timeout,
        )
        if run_chembl(cfg, chembl_args) != 0:
            return 1
        chembl_df = pd.read_csv(
            chembl_out, sep=args.sep, encoding=args.encoding, dtype=str
        ).rename(columns={"target_chembl_id": "chembl_id"})

        # Extract UniProt IDs and write temporary CSV for downstream steps
        uids = [
            u
            for u in chembl_df.get(args.uniprot_column, [])
            if isinstance(u, str) and u
        ]
        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=args.encoding, newline=""
        ) as tmp:
            writer = csv.DictWriter(tmp, fieldnames=["uniprot_id"], delimiter=args.sep)
            writer.writeheader()
            for uid in uids:
                writer.writerow({"uniprot_id": uid})
            tmp_path = Path(tmp.name)

        # Run UniProt pipeline
        uniprot_args = argparse.Namespace(
            input_csv=tmp_path,
            output_csv=uniprot_out,
            data_dir=args.data_dir,
            sep=args.sep,
            encoding=args.encoding,
            column="uniprot_id",
        )
        try:
            if run_uniprot(cfg, uniprot_args) != 0:
                return 1
        finally:
            tmp_path.unlink(missing_ok=True)

        # Load UniProt output
        uniprot_df = pd.read_csv(
            uniprot_out, sep=args.sep, encoding=args.encoding, dtype=str
        )
        # The uids list holds the original identifiers used to query UniProt,
        # and uniprot_df contains the corresponding results in the same order.
        # We add the original IDs to uniprot_df to allow merging with chembl_df.
        uniprot_df["original_id"] = uids

        # To avoid column name collisions during the merge, we drop 'uniprot_id'
        # from chembl_df and rely on the canonical 'uniprot_id' from uniprot_df.
        chembl_for_merge = chembl_df.drop(columns=["uniprot_id"], errors="ignore")

        # Prepare combined input for IUPHAR containing ChEMBL and UniProt data
        combined_df = pd.merge(
            chembl_for_merge,
            uniprot_df,
            left_on=args.uniprot_column,
            right_on="original_id",
            how="left",
        ).drop(columns=["original_id"])

        # Consolidate synonym and EC number information for classification
        combined_df["synonyms"] = combined_df.apply(
            lambda r: _pipe_merge(
                [
                    r.get("pref_name"),
                    r.get("component_description"),
                    r.get("gene"),
                    r.get("chembl_alternative_name"),
                    r.get("names"),
                    r.get("secondaryAccessionNames"),
                ]
            ),
            axis=1,
        )
        combined_df["ec_number"] = combined_df.apply(
            lambda r: _pipe_merge([r.get("ec_numbers"), r.get("reaction_ec_numbers")]),
            axis=1,
        )
        combined_df["gene_name"] = combined_df["gene"].apply(_first_token)
        combined_df = combined_df.drop(
            columns=["ec_numbers", "reaction_ec_numbers"], errors="ignore"
        )

        with NamedTemporaryFile(
            "w", delete=False, encoding=args.encoding, newline=""
        ) as tmp:
            combined_df.to_csv(tmp, index=False, sep=args.sep, encoding=args.encoding)
            iuphar_input = Path(tmp.name)

        # Run IUPHAR mapping using combined data
        iuphar_args = argparse.Namespace(
            input_csv=iuphar_input,
            output_csv=iuphar_out,
            target_csv=args.target_csv,
            family_csv=args.family_csv,
            sep=args.sep,
            encoding=args.encoding,
        )
        try:
            if run_iuphar(cfg, iuphar_args) != 0:
                return 1
        finally:
            iuphar_input.unlink(missing_ok=True)

        # Merge results using pandas

        iuphar_df = pd.read_csv(
            iuphar_out, sep=args.sep, encoding=args.encoding, dtype=str
        )
        existing_cols = set(chembl_df.columns) | set(uniprot_df.columns)
        classification_cols = [c for c in iuphar_df.columns if c not in existing_cols]

        iuphar_df = iuphar_df[["uniprot_id", *classification_cols]]

        merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
        # Apply domain-specific clean-up and finalise table before exporting
        processed = tp.postprocess_targets(merged)
        organism_df = pd.read_csv(
            args.organism_csv, sep=args.sep, encoding=args.encoding, dtype=str
        )
        final_df = tp.finalise_targets(processed, organism_df)
        final_df = normalize_targets(final_df)
        exit_code = 0
        try:
            final_df = TargetsSchema.validate(final_df, lazy=True)
        except SchemaErrors as exc:
            failure_path = output.with_name(f"{output.stem}_failure_cases.csv")
            errors = SidecarErrors()
            for row in exc.failure_cases.to_dict("records"):
                errors.add_error(row)
            errors.save(failure_path)
            logger.error(
                "validation failed; wrote %d failure cases to %s",
                len(exc.failure_cases),
                failure_path,
            )
            final_df = getattr(exc, "validated_data", final_df)
            exit_code = 1
        io.write_csv(final_df, output, cfg=cfg, sep=args.sep, encoding=args.encoding)
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(final_df, table_name=str(output.with_suffix("")))
    except ValueError as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return exit_code


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point using :class:`Config` for defaults."""
    parser, log_cfg = build_parser()
    args = parser.parse_args(argv)
    log_cfg.level = args.log_level
    logger = configure_logger(log_cfg)
    logger.info("pipeline start run_id=%s", log_cfg.run_id, extra={"event": "start"})
    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    try:
        cfg: Config = apply_config_overrides(
            args,
            subparser,
            args.config,
            mapping={
                "timeout": "api.timeout_read",
                "target_csv": "resources.iuphar_target_csv",
                "family_csv": "resources.iuphar_family_csv",
                "data_dir": "resources.uniprot_data_dir",
                "organism_csv": "resources.organism_csv",
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
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
        return 1
    if hasattr(args, "func"):
        exit_code = args.func(cfg, args)
        if exit_code == 0:
            logger.info(
                "pipeline done run_id=%s", log_cfg.run_id, extra={"event": "done"}
            )
        else:
            logger.info(
                "pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"}
            )
        return exit_code
    parser.print_help()
    logger.info("pipeline fail run_id=%s", log_cfg.run_id, extra={"event": "fail"})
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
