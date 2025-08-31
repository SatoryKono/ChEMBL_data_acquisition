"""CLI entry point for review classification."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import List

from review_classifier.pipeline import process_records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify publications as reviews or not")
    parser.add_argument("--input", type=Path, required=True, help="Path to NDJSON input file")
    parser.add_argument("--output", type=Path, required=True, help="Directory for artifacts")
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging level"
    )
    return parser.parse_args()


def load_records(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8") as fh:
        return [json.loads(line) for line in fh]


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper()))
    logger = logging.getLogger("review_classifier")
    records = load_records(args.input)
    process_records(records, args.output, logger)


if __name__ == "__main__":
    main()
