"""Command-line interface for deterministic document classifier."""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime
from typing import Iterable, List

from classifier import DocumentClassifier
from io_schemas import ClassificationResult, DocumentRecord

logger = logging.getLogger(__name__)

LIST_FIELDS = [
    "pubmed_publicationtype",
    "scholar_publicationtypes",
    "openalex_publicationtypes",
    "pubmed_mesh_descriptors",
    "pubmed_mesh_qualifiers",
    "openalex_meshdescriptors",
    "openalex_meshqualifiers",
    "pubmed_chemicallist",
]


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deterministic document classifier")
    parser.add_argument("--input", "--in", dest="input", default="input.csv")
    parser.add_argument("--output", "--out", dest="output")
    parser.add_argument("--in-format", choices=["csv", "jsonl"], default="jsonl")
    parser.add_argument("--out-format", choices=["csv", "jsonl"], default="jsonl")
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def _norm_list_field(value: str) -> List[str]:
    if not value:
        return []
    return [v.strip() for v in value.split("|") if v.strip()]


def read_records(path: str, fmt: str) -> Iterable[DocumentRecord]:
    with open(path, "r", encoding="utf-8") as fh:
        if fmt == "jsonl":
            for line in fh:
                data = json.loads(line)
                for field in LIST_FIELDS:
                    data[field] = _norm_list_field(data.get(field, ""))
                yield DocumentRecord(**data)
        else:
            import csv

            reader = csv.DictReader(fh)
            for row in reader:
                data = {k: v for k, v in row.items()}
                for field in LIST_FIELDS:
                    data[field] = _norm_list_field(data.get(field, ""))
                yield DocumentRecord(**data)


def write_records(records: Iterable[ClassificationResult], path: str, fmt: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        if fmt == "jsonl":
            for rec in records:
                fh.write(rec.model_dump_json() + "\n")
        else:
            import csv

            writer = None
            for rec in records:
                data = rec.model_dump()
                if writer is None:
                    writer = csv.DictWriter(fh, fieldnames=list(data.keys()))
                    writer.writeheader()
                for field in ["vivo_hits", "vitro_hits", "conflicts"]:
                    data[field] = "|".join(data[field])
                assert writer is not None
                writer.writerow(data)


def main() -> None:
    args = _parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    if not args.output:
        dt = datetime.now().strftime("%Y%m%d")
        args.output = f"output_{args.input}_{dt}.{args.out_format}"

    classifier = DocumentClassifier()
    results = [
        classifier.classify(rec) for rec in read_records(args.input, args.in_format)
    ]
    write_records(results, args.output, args.out_format)


if __name__ == "__main__":
    main()
