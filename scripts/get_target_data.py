"""Command line interface for retrieving target data from external sources.

Example
-------
Fetch ChEMBL target information for identifiers in ``targets.csv``::

    python scripts/get_target_data.py chembl --config config.yaml --input targets.csv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import cast

# Allow running the script directly via ``python scripts/get_target_data.py``
# by ensuring the repository root is on ``sys.path`` when the module is executed
# outside of the ``scripts`` package. This mirrors the behaviour of installing
# the project in editable mode.
if __package__ in {None, ""}:
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import pandas as pd
import requests
from pandera.errors import SchemaErrors

from library import chembl_library as cl
from library import io
from library import iuphar_library as ii
from library import target_postprocessing as tp
from library import uniprot_library as uu
from library.chembl_client import ChemblClient
from library.cli import (
    LoggerConfig,
    apply_config_overrides,
    build_root_parser,
    configure_logger,
)
from library.config import (
    Config,
    _serialize_paths,
    ensure_dirs,
    print_config,
)
from library.log import logger
from library.metadata import Stats, file_sha256, write_meta_yaml
from library.sidecar import SidecarErrors
from library.table_quality import analyze_table_quality
from schemas import TargetsSchema, normalize_targets


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


def _save_snapshot(df: pd.DataFrame, base: Path, step: str, cfg: Config) -> Path:
    """Write ``df`` to a uniquely named snapshot CSV file.

    The file is created alongside ``base`` using the pattern
    ``<base>_<step>_<n>.csv`` where ``n`` increments to avoid overwriting
    existing files.

    Parameters
    ----------
    df:
        Data frame to serialise.
    base:
        Base path for the output file. Its stem and suffix determine the
        snapshot file name.
    step:
        Descriptive label inserted into the snapshot file name.
    cfg:
        Application configuration (currently unused but kept for API
        compatibility).

    Returns
    -------
    Path
        Path to the written snapshot file.
    """
    stem = base.stem
    suffix = base.suffix or ".csv"
    index = 1
    while True:
        candidate = base.with_name(f"{stem}_{step}_{index}{suffix}")
        if not candidate.exists():
            df.to_csv(candidate, index=False)
            return candidate
        index += 1


def build_parser() -> tuple[argparse.ArgumentParser, LoggerConfig]:
    """Create and return the top-level CLI argument parser.

    The command line interface is organised into sub-commands for retrieving
    data from individual sources (UniProt, ChEMBL and IUPHAR) as well as a
    convenience ``all`` command that runs all pipelines and merges their
    outputs.
    """

    root, shared, log_cfg = build_root_parser()
    parser = argparse.ArgumentParser(
        description="Target data utilities", parents=[root]
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ----------------------------
    # UniProt sub-command
    # ----------------------------
    uniprot = subparsers.add_parser(
        "uniprot",
        parents=[shared],
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
        parents=[shared],
        help="Retrieve target information from ChEMBL",
    )
    chembl.add_argument(
        "--column",
        default="target_chembl_id",
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
        parents=[shared],
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
        parents=[shared],
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

    parser.subparsers_map = {  # type: ignore[attr-defined]
        "uniprot": uniprot,
        "chembl": chembl,
        "iuphar": iuphar,
        "all": all_cmd,
    }

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
            args.input_csv, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        column = cfg.target.uniprot.column
        if column not in df.columns:
            raise ValueError(f"column '{column}' not found in {args.input_csv}")
        df = df.fillna("")
        df = df[(df[column].str.strip() != "") & (df[column] != "#N/A")].reset_index(
            drop=True
        )
        ids = df[column].tolist()

        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
        ) as tmp:
            writer = csv.DictWriter(
                tmp, fieldnames=["uniprot_id"], delimiter=cfg.io.csv_sep
            )
            writer.writeheader()
            for uid in ids:
                writer.writerow({"uniprot_id": uid})
            tmp_path = Path(tmp.name)

        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        data_dir = cfg.target.uniprot.data_dir
        try:
            uu.process(
                input_csv=str(tmp_path),
                output_csv=str(output),
                data_dir=data_dir,
                cfg=cfg.uniprot,
                sep=cfg.io.csv_sep,
                encoding=cfg.io.csv_encoding,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        out_df = pd.read_csv(
            output, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(1, "mapping_uniprot_id", df["mapping_uniprot_id"].tolist())
        io.write_csv(
            out_df,
            output,
            cfg=cfg,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
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
    # Set up HTTP session with proper headers and retry behaviour
    with ChemblClient(cfg.api, cfg.retry, cfg.chembl) as client:
        try:
            ids = io.read_ids(
                args.input_csv, column=cfg.target.chembl.column, cfg=cfg.io
            )
        except (FileNotFoundError, ValueError) as exc:
            logger.error("%s", exc)
            return 1

        try:
            df = cl.get_targets(
                ids,
                cfg=cfg.api,
                client=client,
                mapping_cfg=cfg.uniprot_mapping,
                timeout=cfg.target.chembl.timeout,
            )
        except (requests.RequestException, ValueError) as exc:
            logger.error("failed to retrieve targets: %s", exc)
            return 1
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        df = normalize_targets(df)
        rows_total = len(df)
        exit_code = 0
    # Only enforce columns marked as ``required`` in the schema. Optional columns
    # may be absent in the input dataset without preventing validation.
    required_cols = {
        name for name, col in TargetsSchema.columns.items() if col.required
    }
    optional_cols = set(TargetsSchema.columns) - required_cols
    missing_required = required_cols - set(df.columns)
    missing_optional = optional_cols - set(df.columns)
    if not missing_required:
        if missing_optional:
            logger.warning(
                "DataFrame is missing optional columns: %s", missing_optional
            )
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
        logger.warning(
            "Skipping validation due to missing required columns: %s",
            missing_required,
        )
    rows_kept = len(df)
    rows_dropped = rows_total - rows_kept
    try:
        key_cols = [c for c in ["target_chembl_id"] if c in df.columns]
        csv_path = io.write_csv(
            df,
            output,
            cfg=cfg,
            key_cols=key_cols or None,
        )
        logger.info("write_done", extra={"rows": rows_kept, "path": str(csv_path)})
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
            target_path=cfg.target.iuphar.target_csv,
            family_path=cfg.target.iuphar.family_csv,
            encoding=cfg.io.csv_encoding,
        )
        output = args.output_csv or io.default_output_path(args.input_csv, cfg.io)
        data.map_uniprot_file(
            input_path=args.input_csv,
            output_path=output,
            encoding=cfg.io.csv_encoding,
            sep=cfg.io.csv_sep,
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
        chembl_out = cfg.target.all.chembl_out or output.with_name(
            output.stem + "_chembl.csv"
        )
        uniprot_out = cfg.target.all.uniprot_out or output.with_name(
            output.stem + "_uniprot.csv"
        )
        iuphar_out = cfg.target.all.iuphar_out or output.with_name(
            output.stem + "_iuphar.csv"
        )

        # Run ChEMBL retrieval and capture results
        chembl_args = argparse.Namespace(
            input_csv=args.input_csv,
            output_csv=chembl_out,
        )
        if run_chembl(cfg, chembl_args) != 0:
            return 1
        chembl_df = pd.read_csv(
            chembl_out, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        ).rename(columns={"target_chembl_id": "target_chembl_id"})

        # Extract UniProt IDs and write temporary CSV for downstream steps
        uids = [
            u
            for u in chembl_df.get(cfg.target.all.uniprot_column, [])
            if isinstance(u, str) and u
        ]
        from tempfile import NamedTemporaryFile

        with NamedTemporaryFile(
            "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
        ) as tmp:
            writer = csv.DictWriter(
                tmp, fieldnames=["uniprot_id"], delimiter=cfg.io.csv_sep
            )
            writer.writeheader()
            for uid in uids:
                writer.writerow({"uniprot_id": uid})
            tmp_path = Path(tmp.name)

        # Run UniProt pipeline
        uniprot_args = argparse.Namespace(
            input_csv=tmp_path,
            output_csv=uniprot_out,
        )
        orig_dir = cfg.target.uniprot.data_dir
        cfg.target.uniprot.data_dir = cfg.target.all.data_dir
        try:
            if run_uniprot(cfg, uniprot_args) != 0:
                return 1
        finally:
            cfg.target.uniprot.data_dir = orig_dir
            tmp_path.unlink(missing_ok=True)

        # Load UniProt output
        uniprot_df = pd.read_csv(
            uniprot_out, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        # The uids list holds the original identifiers used to query UniProt,
        # and uniprot_df contains the corresponding results in the same order.
        # We add the original IDs to uniprot_df to allow merging with chembl_df.
        uniprot_df["original_id"] = uids

        # To avoid column name collisions during the merge, drop "uniprot_id"
        # from the ChEMBL frame only when it is *not* used as the merge key.
        # When "uniprot_id" is the join column we keep it for the merge and
        # replace it afterwards with the canonical value from ``uniprot_df``.
        if cfg.target.all.uniprot_column != "uniprot_id":
            chembl_for_merge = chembl_df.drop(columns=["uniprot_id"], errors="ignore")
        else:
            chembl_for_merge = chembl_df.copy()

        # Prepare combined input for IUPHAR containing ChEMBL and UniProt data.
        combined_df = pd.merge(
            chembl_for_merge,
            uniprot_df,
            left_on=cfg.target.all.uniprot_column,
            right_on="original_id",
            how="left",
        ).drop(columns=["original_id"])
        if cfg.target.all.uniprot_column == "uniprot_id":
            # ``uniprot_df`` contains the authoritative ID, so discard the one
            # from ChEMBL and rename the UniProt column back to ``uniprot_id``.
            combined_df = combined_df.drop(columns=["uniprot_id_x"], errors="ignore")
            combined_df = combined_df.rename(columns={"uniprot_id_y": "uniprot_id"})

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
            "w", delete=False, encoding=cfg.io.csv_encoding, newline=""
        ) as tmp:
            combined_df.to_csv(
                tmp, index=False, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding
            )
            iuphar_input = Path(tmp.name)

        # Run IUPHAR mapping using combined data
        iuphar_args = argparse.Namespace(
            input_csv=iuphar_input,
            output_csv=iuphar_out,
        )
        orig_target = cfg.target.iuphar.target_csv
        orig_family = cfg.target.iuphar.family_csv
        cfg.target.iuphar.target_csv = cfg.target.all.target_csv
        cfg.target.iuphar.family_csv = cfg.target.all.family_csv
        try:
            if run_iuphar(cfg, iuphar_args) != 0:
                return 1
        finally:
            cfg.target.iuphar.target_csv = orig_target
            cfg.target.iuphar.family_csv = orig_family
            iuphar_input.unlink(missing_ok=True)

        # Merge results using pandas

        iuphar_df = pd.read_csv(
            iuphar_out, sep=cfg.io.csv_sep, encoding=cfg.io.csv_encoding, dtype=str
        )
        # ``combined_df`` already contains the consolidated data from
        # ChEMBL and UniProt, including helper fields like ``synonyms``,
        # ``ec_number`` and ``gene_name``. Build the set of existing
        # columns from this frame to avoid reintroducing those fields
        # when selecting IUPHAR classification columns.
        existing_cols = set(combined_df.columns)
        classification_cols = [c for c in iuphar_df.columns if c not in existing_cols]

        # Recreate the IUPHAR frame with only the new classification
        # columns before merging.
        iuphar_df = iuphar_df[["uniprot_id", *classification_cols]].copy()

        merged = combined_df.merge(iuphar_df, on="uniprot_id", how="left")
        # Apply domain-specific clean-up and finalise table before exporting
        processed = tp.postprocess_targets(merged)
        organism_df = pd.read_csv(
            cfg.target.all.organism_csv,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
            dtype=str,
        )
        final_df = tp.finalise_targets(processed, organism_df)
        # ``postprocess_targets`` and ``finalise_targets`` operate on the
        # ``chembl_id`` column name, but the validation schema expects the
        # original ``target_chembl_id``. Rename the column back before
        # normalisation and validation to satisfy the schema requirements.
       # final_df = final_df.rename(columns={"chembl_id": "target_chembl_id"})
        final_df = normalize_targets(final_df)
        exit_code = 0
        required_cols = {
            name for name, col in TargetsSchema.columns.items() if col.required
        }
        optional_cols = set(TargetsSchema.columns) - required_cols
        missing_required = required_cols - set(final_df.columns)
        missing_optional = optional_cols - set(final_df.columns)
        if not missing_required:
            if missing_optional:
                logger.warning(
                    "DataFrame is missing optional columns: %s", missing_optional
                )
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
        else:
            logger.warning(
                "Skipping validation due to missing required columns: %s",
                missing_required,
            )
        final_df = final_df.drop_duplicates()
        io.write_csv(
            final_df,
            output,
            cfg=cfg,
            sep=cfg.io.csv_sep,
            encoding=cfg.io.csv_encoding,
        )
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
    logger.info("pipeline_start", extra={"run_id": log_cfg.run_id})
    subparser_map = getattr(parser, "subparsers_map", {})
    subparser = subparser_map.get(args.command, parser)
    try:
        mapping: dict[str, str] = {}
        if args.command == "uniprot":
            mapping = {
                "column": "target.uniprot.column",
                "data_dir": "target.uniprot.data_dir",
            }
        elif args.command == "chembl":
            mapping = {
                "column": "target.chembl.column",
                "timeout": "target.chembl.timeout",
            }
        elif args.command == "iuphar":
            mapping = {
                "target_csv": "target.iuphar.target_csv",
                "family_csv": "target.iuphar.family_csv",
            }
        elif args.command == "all":
            mapping = {
                "timeout": "target.all.timeout",
                "data_dir": "target.all.data_dir",
                "target_csv": "target.all.target_csv",
                "family_csv": "target.all.family_csv",
                "organism_csv": "target.all.organism_csv",
                "uniprot_column": "target.all.uniprot_column",
                "chembl_out": "target.all.chembl_out",
                "uniprot_out": "target.all.uniprot_out",
                "iuphar_out": "target.all.iuphar_out",
            }
        cfg: Config = apply_config_overrides(
            args, subparser, args.config, mapping=mapping, base_parser=parser
        )
        if args.print_config:
            print_config(cfg)
            configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
            logger.info("pipeline_done", extra={"run_id": log_cfg.run_id})
            return 0
        ensure_dirs(cfg)
        logger = configure_logger(log_cfg, fmt=cfg.log.format, datefmt=cfg.log.datefmt)
    except (ValueError, TypeError) as exc:
        logger.error("%s", exc)
        logger.info("pipeline_fail", extra={"run_id": log_cfg.run_id})
        return 1
    except (FileNotFoundError, NotADirectoryError) as exc:
        logger.error("failed to set up directories: %s", exc)
        logger.info("pipeline_fail", extra={"run_id": log_cfg.run_id})
        return 1
    if hasattr(args, "func"):
        exit_code = cast(int, args.func(cfg, args))
        if exit_code == 0:
            logger.info("pipeline_done", extra={"run_id": log_cfg.run_id})
        else:
            logger.info("pipeline_fail", extra={"run_id": log_cfg.run_id})
        return exit_code
    parser.print_help()
    logger.info("pipeline_fail", extra={"run_id": log_cfg.run_id})
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
