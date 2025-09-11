"""CLI wrapper for :func:`write_csv_deterministic`.

This script reads an input CSV file and re-serialises it deterministically using
:func:`chembl_da.library.csv_utils.write_csv_deterministic`.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import pandas as pd

from chembl_da.library.csv_utils import write_csv_deterministic


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="input.csv", help="Input CSV path")
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (default: derive from input name)",
    )
    parser.add_argument("--col-order", nargs="*", help="Preferred column order")
    parser.add_argument("--key-cols", nargs="*", help="Columns used for sorting")
    parser.add_argument("--sep", default=",", help="Input CSV delimiter")
    parser.add_argument(
        "--encoding",
        default="utf8",
        help="Input CSV encoding",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    df = pd.read_csv(args.input, sep=args.sep, encoding=args.encoding)
    output = args.output or Path(args.input).with_name(
        f"output_{Path(args.input).stem}.csv"
    )
    write_csv_deterministic(
        df,
        output,
        col_order=args.col_order or None,
        key_cols=args.key_cols or None,
    )
    logging.info("written %s", output)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
