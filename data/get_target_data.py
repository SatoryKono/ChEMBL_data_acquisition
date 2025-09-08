"""Command line interface for retrieving target data from external sources."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import chembl_library as cl
from library import io
from library import iuphar_library as ii
from library import target_postprocessing as tp
from library import uniprot_library as uu
from library import utils
from library.table_quality import analyze_table_quality

logger = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    """Create and return the top-level CLI argument parser.

    The command line interface is organised into sub-commands for
    retrieving data from individual sources (UniProt, ChEMBL and IUPHAR)
    as well as a convenience ``all`` command that runs all pipelines and
    merges their outputs.
    """

    parser = argparse.ArgumentParser(description="Target data utilities")
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING)",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # ----------------------------
    # UniProt sub-command
    # ----------------------------
    uniprot = subparsers.add_parser(
        "uniprot", help="Extract information for UniProt accessions"
    )
    uniprot.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV file containing a 'uniprot_id' column",
    )
    uniprot.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file for the extracted information (default: auto-generate)",
    )
    uniprot.add_argument(
        "--column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column in the input CSV containing UniProt accessions",
    )
    uniprot.add_argument(
        "--sep",
        default=",",
        help="CSV delimiter used for input and output files",
    )
    uniprot.add_argument(
        "--encoding",
        default="utf8",
        help="File encoding for input and output CSV files",
    )
    uniprot.add_argument(
        "--data-dir",
        default="uniprot",
        help="Directory containing '<uniprot_id>.json' files",
    )
    uniprot.set_defaults(func=run_uniprot)

    # ----------------------------
    # ChEMBL sub-command
    # ----------------------------
    chembl = subparsers.add_parser(
        "chembl", help="Retrieve target information from ChEMBL"
    )
    chembl.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV file containing ChEMBL target identifiers",
    )
    chembl.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file for target information (default: auto-generate)",
    )
    chembl.add_argument(
        "--column",
        default="chembl_id",
        help="Column name in the input CSV containing identifiers",
    )
    chembl.add_argument("--sep", default=",", help="CSV delimiter for I/O")
    chembl.add_argument(
        "--encoding",
        default="utf8",
        help="File encoding for input and output CSV files",
    )
    chembl.set_defaults(func=run_chembl)

    # ----------------------------
    # IUPHAR sub-command
    # ----------------------------
    iuphar = subparsers.add_parser(
        "iuphar", help="Map UniProt accessions to IUPHAR classifications"
    )
    iuphar.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV file containing a 'uniprot_id' column",
    )
    iuphar.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file for the mapping results (default: auto-generate)",
    )
    iuphar.add_argument(
        "--target-csv",
        type=Path,
        default=Path("dictionary/_IUPHAR/_IUPHAR_target.csv"),
        help="Path to the _IUPHAR_target.csv file",
    )
    iuphar.add_argument(
        "--family-csv",
        type=Path,
        default=Path("dictionary/_IUPHAR/_IUPHAR_family.csv"),
        help="Path to the _IUPHAR_family.csv file",
    )
    iuphar.add_argument("--sep", default=",", help="CSV delimiter for I/O")
    iuphar.add_argument(
        "--encoding",
        default="utf8",
        help="File encoding for input and output CSV files",
    )
    iuphar.set_defaults(func=run_iuphar)

    # ----------------------------
    # Combined pipeline
    # ----------------------------
    all_cmd = subparsers.add_parser(
        "all",
        help="Run ChEMBL, UniProt and IUPHAR pipelines and merge results",
    )
    all_cmd.add_argument(
        "--input",
        dest="input_csv",
        type=Path,
        default=Path("input.csv"),
        help="CSV with a 'chembl_id' column",
    )
    all_cmd.add_argument(
        "--output",
        dest="output_csv",
        type=Path,
        default=None,
        help="Destination CSV file for the merged table (default: auto-generate)",
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
        default="uniprot",
        help="Directory containing '<uniprot_id>.json' files",
    )
    all_cmd.add_argument(
        "--target-csv",
        type=Path,
        default=Path("dictionary/_IUPHAR/_IUPHAR_target.csv"),
        help="Path to the _IUPHAR_target.csv file",
    )
    all_cmd.add_argument(
        "--family-csv",
        type=Path,
        default=Path("dictionary/_IUPHAR/_IUPHAR_family.csv"),
        help="Path to the _IUPHAR_family.csv file",
    )
    all_cmd.add_argument(
        "--organism-csv",
        type=Path,
        default=Path("dictionary/organism.csv"),
        help="CSV mapping 'genus' to organism 'type' for finalisation",
    )
    all_cmd.add_argument(
        "--uniprot-column",
        default="uniprot_id",
        choices=["uniprot_id", "mapping_uniprot_id"],
        help="Column from ChEMBL output to use for UniProt processing",
    )
    all_cmd.add_argument("--sep", default=",", help="CSV delimiter for I/O")
    all_cmd.add_argument(
        "--encoding",
        default="utf8",
        help="File encoding for input and output CSV files",
    )
    all_cmd.set_defaults(func=run_all)

    return parser


def configure_logging(level: str) -> None:
    """Configure basic logging.

    Parameters
    ----------
    level:
        Logging verbosity such as ``"INFO"`` or ``"DEBUG"``.

    Returns
    -------
    None
        This function configures the root logger in-place and does not
        return a value.
    """
    logging.basicConfig(level=getattr(logging, level.upper(), logging.INFO))


def run_uniprot(args: argparse.Namespace) -> int:
    """Execute the ``uniprot`` sub-command.

    Parameters
    ----------
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

        output = args.output_csv or io.default_output_path(args.input_csv)
        uu.process_ids(
            uniprot_ids=ids,
            output_csv=str(output),
            data_dir=str(args.data_dir),
            sep=args.sep,
            encoding=args.encoding,
        )

        out_df = pd.read_csv(output, sep=args.sep, encoding=args.encoding, dtype=str)
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(1, "mapping_uniprot_id", df["mapping_uniprot_id"].tolist())
        io.write_csv(out_df, output, sep=args.sep, encoding=args.encoding)
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    try:
        analyze_table_quality(out_df, table_name=str(output.with_suffix("")))
    except Exception as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_chembl(args: argparse.Namespace) -> int:
    """Execute the ``chembl`` sub-command."""

    try:
        ids = io.read_ids(
            args.input_csv, column=args.column, sep=args.sep, encoding=args.encoding
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    df = cl.get_targets(ids)
    output = args.output_csv or io.default_output_path(args.input_csv)
    try:
        io.write_csv(df, output, sep=args.sep, encoding=args.encoding)
        logger.info("Wrote %d rows to %s", len(df), output)
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
        return 1
    try:
        analyze_table_quality(df, table_name=str(output.with_suffix("")))
    except Exception as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def run_iuphar(args: argparse.Namespace) -> int:
    """Execute the ``iuphar`` sub-command."""

    try:
        data = ii.IUPHARData.from_files(
            target_path=args.target_csv,
            family_path=args.family_csv,
            encoding=args.encoding,
        )
        output = args.output_csv or io.default_output_path(args.input_csv)
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
    except Exception as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def _fetch_chembl_data(args: argparse.Namespace) -> pd.DataFrame:
    """Fetch ChEMBL data."""
    ids = io.read_ids(
        args.input_csv, column="chembl_id", sep=args.sep, encoding=args.encoding
    )
    df = cl.get_targets(ids)
    return df.rename(columns={"target_chembl_id": "chembl_id"})


def _fetch_uniprot_data(
    uids: list[str], args: argparse.Namespace
) -> pd.DataFrame:
    """Fetch UniProt data."""
    return uu.get_info_for_ids(uids, data_dir=args.data_dir)


def _merge_chembl_uniprot(
    chembl_df: pd.DataFrame, uniprot_df: pd.DataFrame, uniprot_column: str
) -> pd.DataFrame:
    """Merge ChEMBL and UniProt data."""
    uniprot_df["original_id"] = uniprot_df["uniprot_id"]
    chembl_for_merge = chembl_df.drop(columns=["uniprot_id"], errors="ignore")
    combined_df = pd.merge(
        chembl_for_merge,
        uniprot_df,
        left_on=uniprot_column,
        right_on="original_id",
        how="left",
    ).drop(columns=["original_id"])
    return combined_df


def _classify_with_iuphar(
    df: pd.DataFrame, args: argparse.Namespace
) -> pd.DataFrame:
    """Classify the data with IUPHAR."""
    iuphar_data = ii.IUPHARData.from_files(
        target_path=args.target_csv,
        family_path=args.family_csv,
        encoding=args.encoding,
    )
    return iuphar_data.map_uniprot_df(df)


def _postprocess_and_finalize(
    df: pd.DataFrame, args: argparse.Namespace
) -> pd.DataFrame:
    """Post-process and finalize the data."""
    df["synonyms"] = df.apply(
        lambda r: utils.pipe_merge(
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
    df["ec_number"] = df.apply(
        lambda r: utils.pipe_merge(
            [r.get("ec_numbers"), r.get("reaction_ec_numbers")]
        ),
        axis=1,
    )
    df["gene_name"] = df["gene"].apply(utils.first_token)
    df = df.drop(columns=["ec_numbers", "reaction_ec_numbers"], errors="ignore")

    processed = tp.postprocess_targets(df)
    organism_df = pd.read_csv(
        args.organism_csv, sep=args.sep, encoding=args.encoding, dtype=str
    )
    return tp.finalise_targets(processed, organism_df)


def run_all(args: argparse.Namespace) -> int:
    """Run ChEMBL, UniProt and IUPHAR pipelines and merge their outputs."""
    try:
        chembl_df = _fetch_chembl_data(args)
        uids = [
            u
            for u in chembl_df.get(args.uniprot_column, [])
            if isinstance(u, str) and u
        ]
        uniprot_df = _fetch_uniprot_data(uids, args)
        combined_df = _merge_chembl_uniprot(
            chembl_df, uniprot_df, args.uniprot_column
        )
        classified_df = _classify_with_iuphar(combined_df, args)
        final_df = _postprocess_and_finalize(classified_df, args)

        output = args.output_csv or io.default_output_path(args.input_csv)
        io.write_csv(final_df, output, sep=args.sep, encoding=args.encoding)

        analyze_table_quality(final_df, table_name=str(output.with_suffix("")))

    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1
    except Exception as exc:
        logger.error("failed to generate quality report: %s", exc)
        return 1
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    """Command line entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    if hasattr(args, "func"):
        return args.func(args)
    parser.print_help()
    return 1


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
