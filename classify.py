"""CLI for document review classification."""
from __future__ import annotations

import argparse
import csv
import json
import logging
from pathlib import Path
from typing import List, Dict, Any

from doc_classifier import io_utils, normalize, signals, vote, logging_utils
from doc_classifier.weights import DEFAULT_THRESHOLD


def process_file(
    input_path: str,
    output_path: str,
    log_path: str | None = None,
    chunk_size: int = 1000,
    threshold: float = DEFAULT_THRESHOLD,
) -> None:
    """Process the input file and write classification results."""

    out_fields = [
        "document_chembl_id",
        "doi",
        "pubmed_id",
        "openalex.paperid",
        "label",
        "review_score",
        "nonreview_score",
        "decision_margin",
        "decision_path",
    ]

    logger = logging.getLogger(__name__)
    logger.info("Starting classification of %s", input_path)

    # Prepare output and log files
    with open(output_path, "w", newline="", encoding="utf-8") as out_fp:
        writer = csv.DictWriter(out_fp, fieldnames=out_fields)
        writer.writeheader()

        log_fp = open(log_path, "w", encoding="utf-8") if log_path else None
        try:
            row_index = -1
            for chunk in io_utils.read_chunks(input_path, chunk_size=chunk_size):
                for _, row in chunk.iterrows():
                    row_index += 1
                    ids = {col: row.get(col, "") for col in io_utils.ID_COLUMNS}

                    normalized: Dict[str, List[str]] = {}
                    fires: List[signals.Signal] = []
                    notes: List[str] = []
                    protocol_seen = False

                    # PublicationType fields
                    for source, col in [
                        ("pubmed", "pubmed_pt"),
                        ("openalex", "openalex_pt"),
                        ("scholar", "scholar_pt"),
                        ("crossref", "crossref_type"),
                    ]:
                        tokens = normalize.split_and_canon(row.get(col))
                        normalized[col] = tokens
                        if "journal-article" in tokens:
                            notes.append("journal-article seen without effect")
                        if "protocol" in tokens:
                            protocol_seen = True
                        fires.extend(signals.extract_pt_signals(tokens, source))

                    # MeSH fields
                    for source, desc_col, qual_col in [
                        ("pubmed", "pubmed_mesh_desc", "pubmed_mesh_qual"),
                        (
                            "openalex",
                            "openalex_mesh_desc",
                            "openalex_mesh_qual",
                        ),
                    ]:
                        desc_tokens = normalize.split_and_canon(row.get(desc_col))
                        qual_tokens = normalize.split_and_canon(row.get(qual_col))
                        normalized[desc_col] = desc_tokens
                        normalized[qual_col] = qual_tokens
                        fires.extend(
                            signals.extract_mesh_signals(
                                desc_tokens, qual_tokens, source
                            )
                        )

                    scores = vote.score(fires)
                    label = vote.decide(scores, threshold)
                    if protocol_seen:
                        label = "unknown"
                        notes.append("protocol → unknown rule applied? true")
                    else:
                        notes.append("protocol → unknown rule applied? false")

                    decision_margin = scores["review"] - scores["nonreview"]

                    def sig_to_dict(sig: signals.Signal | None) -> Dict[str, Any] | None:
                        if sig is None:
                            return None
                        return {
                            "source": sig.source,
                            "field": sig.field,
                            "kind": sig.kind,
                            "token": sig.token,
                            "points": sig.points,
                        }

                    top_review = max(
                        (s for s in fires if s.kind == "review"),
                        key=lambda s: s.points,
                        default=None,
                    )
                    top_nonreview = max(
                        (s for s in fires if s.kind == "non-review"),
                        key=lambda s: s.points,
                        default=None,
                    )
                    decision_path = {
                        "top_review": sig_to_dict(top_review),
                        "top_nonreview": sig_to_dict(top_nonreview),
                    }

                    out_row = {
                        **ids,
                        "label": label,
                        "review_score": round(scores["review"], 3),
                        "nonreview_score": round(scores["nonreview"], 3),
                        "decision_margin": round(decision_margin, 3),
                        "decision_path": json.dumps(
                            decision_path, ensure_ascii=False, separators=(",", ":")
                        ),
                    }
                    writer.writerow(out_row)

                    if log_fp:
                        log_record = {
                            "row_index": row_index,
                            "ids": ids,
                            "normalized": normalized,
                            "fires": [sig.__dict__ for sig in fires],
                            "scores": scores,
                            "label": label,
                            "notes": notes,
                        }
                        logging_utils.emit_log(log_record, log_fp)
        finally:
            if log_fp:
                log_fp.close()
    logger.info("Finished. Results written to %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Document review classifier")
    parser.add_argument("--input", required=True, help="Input CSV/TSV file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument("--log", help="Optional JSONL log file")
    parser.add_argument("--chunk-size", type=int, default=1000)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    process_file(
        input_path=args.input,
        output_path=args.output,
        log_path=args.log,
        chunk_size=args.chunk_size,
        threshold=args.threshold,
    )


if __name__ == "__main__":
    main()
