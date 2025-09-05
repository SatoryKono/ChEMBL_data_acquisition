"""Command line interface for retrieving target data from external sources."""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from library import chembl_library as cl
from library import iuphar_library as ii
from library import uniprot_library as uu
from library import io

logger = logging.getLogger(__name__)


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

        output = args.output_csv or io.default_output_path(args.input_csv)
        try:
            uu.process(
                input_csv=str(tmp_path),
                output_csv=str(output),
                data_dir=str(args.data_dir),
                sep=args.sep,
                encoding=args.encoding,
            )
        finally:
            tmp_path.unlink(missing_ok=True)

        out_df = pd.read_csv(output, sep=args.sep, encoding=args.encoding, dtype=str)
        if "mapping_uniprot_id" in df.columns:
            out_df.insert(1, "mapping_uniprot_id", df["mapping_uniprot_id"].tolist())
        io.write_csv(out_df, output, sep=args.sep, encoding=args.encoding)
        return 0
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1


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
    except OSError as exc:
        logger.error("failed to write output CSV: %s", exc)
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
        return 0
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1


def run_all(args: argparse.Namespace) -> int:
    """Run ChEMBL, UniProt and IUPHAR pipelines and merge their outputs.

    Parameters
    ----------
    args:
        Parsed command-line arguments specific to the ``all`` sub-command.

    Returns
    -------
    int
        Zero on success, non-zero on failure.
    """

    try:
        output = args.output_csv or io.default_output_path(args.input_csv)
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
        )
        if run_chembl(chembl_args) != 0:
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
            if run_uniprot(uniprot_args) != 0:
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
            if run_iuphar(iuphar_args) != 0:
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

        io.write_csv(merged, output, sep=args.sep, encoding=args.encoding)
        return 0
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1


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
