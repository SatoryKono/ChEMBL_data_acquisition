"""CLI entry point for review classification."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import List

from review_classifier.pipeline import process_records


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Classify publications as reviews or not")
    parser.add_argument("--input", type=Path, required=True, help="Path to NDJSON input file")
    parser.add_argument("--output", type=Path, required=True, help="Directory for artifacts")
    parser.add_argument("--encoding", default="utf-8", help="Encoding of the input file")
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level"
    )
    return parser.parse_args()


def load_records(path: Path, encoding: str) -> List[dict]:
    """Load NDJSON records from ``path`` using ``encoding``.

    Parameters
    ----------
    path : Path
        Location of the NDJSON file.
    encoding : str
        Text encoding used by ``path``.

    Returns
    -------
    list of dict
        Parsed records from the input file.

    Raises
    ------
    ValueError
        If the file cannot be decoded or contains malformed JSON.
    FileNotFoundError
        If ``path`` does not exist.
    """
    try:
        with path.open("r", encoding=encoding) as fh:
            return [json.loads(line) for line in fh]
    except UnicodeDecodeError as exc:
        raise ValueError(f"failed to decode {path} with encoding '{encoding}'") from exc
    except json.JSONDecodeError as exc:
        raise ValueError(f"malformed JSON in {path}: {exc}") from exc
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"input file not found: {path}") from exc


def main() -> None:
    """Command-line interface for document classification."""
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    logger = logging.getLogger("review_classifier")
    try:
        records = load_records(args.input, args.encoding)
    except (ValueError, FileNotFoundError) as exc:
        logger.error("%s", exc)
        raise SystemExit(1)
    process_records(records, args.output, logger)


if __name__ == "__main__":
    main()
