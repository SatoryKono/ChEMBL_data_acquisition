"""CLI for document type classification."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from mylib.io_utils import ParquetBatchWriter, read_input
from mylib import transforms
from mylib.validate import ensure_columns


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify documents as review or non-review")
    parser.add_argument("--input", required=True, type=Path, help="Path to input CSV/TSV file")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory to store artifacts")
    parser.add_argument("--sep", default=",", help="Field separator of the input file")
    parser.add_argument("--encoding", default="utf-8", help="File encoding")
    parser.add_argument("--batch-size", type=int, default=1000, help="Number of rows per batch")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))

    df = read_input(args.input, args.sep, args.encoding)
    ensure_columns(df)

    out_dir = args.output_dir
    (out_dir / "normalized").mkdir(parents=True, exist_ok=True)
    (out_dir / "features").mkdir(parents=True, exist_ok=True)
    (out_dir / "scores").mkdir(parents=True, exist_ok=True)
    (out_dir / "decisions").mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)

    pt_writer = ParquetBatchWriter(out_dir / "normalized" / "publication_type.parquet")
    mesh_writer = ParquetBatchWriter(out_dir / "features" / "mesh.parquet")
    score_writer = ParquetBatchWriter(out_dir / "scores" / "scores.parquet")
    label_writer = ParquetBatchWriter(out_dir / "decisions" / "labels.parquet")
    log_path = out_dir / "logs" / "decision.log"
    with log_path.open("w", encoding="utf-8") as log_file:
        total = len(df)
        for start in range(0, total, args.batch_size):
            batch = df.iloc[start : start + args.batch_size]
            pt_rows = []
            mesh_rows = []
            score_rows = []
            label_rows = []
            for _, row in batch.iterrows():
                res = transforms.classify_row(row, log_file)
                pt_rows.append({
                    "record_id": res.record_id,
                    **{f"{k}_pt": ";".join(v) for k, v in res.normalized_pt.items()},
                })
                mesh_rows.append({
                    "record_id": res.record_id,
                    **{k: ";".join(v) for k, v in res.normalized_mesh.items()},
                })
                score_rows.append({
                    "record_id": res.record_id,
                    "score_review": res.score_review,
                    "score_non_review": res.score_non_review,
                })
                label_rows.append({
                    "record_id": res.record_id,
                    "label": res.label,
                    "reason": res.reason,
                    "score_review": res.score_review,
                    "score_non_review": res.score_non_review,
                    "factors": ";".join(res.factors),
                })
            pt_writer.write(pd.DataFrame(pt_rows))
            mesh_writer.write(pd.DataFrame(mesh_rows))
            score_writer.write(pd.DataFrame(score_rows))
            label_writer.write(pd.DataFrame(label_rows))
            logging.info("Processed %s/%s records", min(start + args.batch_size, total), total)
    pt_writer.close()
    mesh_writer.close()
    score_writer.close()
    label_writer.close()

    meta = {"rules_version": transforms.RE_RULES_VERSION, "pt_mapping_version": "1"}
    with (out_dir / "run_meta.json").open("w", encoding="utf-8") as fh:
        json.dump(meta, fh)


if __name__ == "__main__":
    main()
