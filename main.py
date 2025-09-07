"""Command-line interface for deterministic document classifier."""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List

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

CSV_FIELD_MAP: Dict[str, str] = {
    "PubMed.PublicationType": "pubmed_publicationtype",
    "scholar.PublicationTypes": "scholar_publicationtypes",
    "OpenAlex.PublicationTypes": "openalex_publicationtypes",
    "crossref.Type": "crossref_type",
    "OpenAlex.TypeCrossref": "openalex_typecrossref",
    "PubMed.MeSH_Descriptors": "pubmed_mesh_descriptors",
    "PubMed.MeSH_Qualifiers": "pubmed_mesh_qualifiers",
    "OpenAlex.MeshDescriptors": "openalex_meshdescriptors",
    "OpenAlex.MeshQualifiers": "openalex_meshqualifiers",
    "PubMed.ChemicalList": "pubmed_chemicallist",
}


def _parse_args() -> argparse.Namespace:
    """Return command-line arguments."""

    parser = argparse.ArgumentParser(description="Deterministic document classifier")
    parser.add_argument("--input", "--in", dest="input", default="input.csv")
    parser.add_argument("--output", "--out", dest="output")
    parser.add_argument("--in-format", choices=["csv", "jsonl"], default=None)
    parser.add_argument("--out-format", choices=["csv", "jsonl"], default=None)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def _norm_list_field(value: str) -> List[str]:
    if not value:
        return []
    return [v.strip() for v in value.split("|") if v.strip()]


def read_records(path: str, fmt: str) -> Iterator[DocumentRecord]:
    """Yield records from ``path`` according to ``fmt``."""

    with open(path, "r", encoding="utf-8") as fh:
        if fmt == "jsonl":
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                except json.JSONDecodeError as exc:  # pragma: no cover - defensive
                    raise ValueError(
                        "Failed to parse JSON; specify --in-format csv for CSV input"
                    ) from exc
                for field in LIST_FIELDS:
                    data[field] = _norm_list_field(data.get(field, ""))
                yield DocumentRecord(**data)
        else:
            import csv

            reader = csv.DictReader(fh)
            for row in reader:
                data = {CSV_FIELD_MAP.get(k, k): v for k, v in row.items()}
                for field in LIST_FIELDS:
                    data[field] = _norm_list_field(data.get(field, ""))
                yield DocumentRecord(**data)


def write_records(records: Iterable[ClassificationResult], path: str, fmt: str) -> None:
    """Write classification ``records`` to ``path`` in ``fmt`` format."""

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
    """CLI entry point."""

    args = _parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    if args.in_format is None:
        args.in_format = "csv" if args.input.lower().endswith(".csv") else "jsonl"

    if args.output:
        if args.out_format is None:
            args.out_format = "csv" if args.output.lower().endswith(".csv") else "jsonl"
    else:
        dt = datetime.now().strftime("%Y%m%d")
        stem = Path(args.input).stem
        out_fmt = args.out_format or args.in_format
        args.output = f"output_{stem}_{dt}.{out_fmt}"
        if args.out_format is None:
            args.out_format = out_fmt

    classifier = DocumentClassifier()
    records = read_records(args.input, args.in_format)
    results = (classifier.classify(rec) for rec in records)
    write_records(results, args.output, args.out_format)


if __name__ == "__main__":
    main()
