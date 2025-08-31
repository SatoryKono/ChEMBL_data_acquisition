from __future__ import annotations

"""CLI for classifying documents as review or non-review."""

import argparse
import logging
from pathlib import Path
from typing import Any, Dict

import yaml

from review_classifier.pipeline import Config, run_pipeline


LOGGER = logging.getLogger("review_classifier")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify documents using publication types and MeSH")
    parser.add_argument("--input-all", dest="input_all", type=Path, required=True, help="Path to all_documents.csv")
    parser.add_argument("--input-mesh", dest="input_mesh", type=Path, required=True, help="Path to experimental_mesh.csv")
    parser.add_argument("--output", dest="output", type=Path, required=True, help="Output file path")
    parser.add_argument("--config", dest="config", type=Path, default=None, help="Optional YAML configuration file")
    parser.add_argument("--delta", dest="delta", type=float, default=0.5, help="Delta threshold for MeSH refinement")
    parser.add_argument("--k-min", dest="k_min", type=int, default=3, help="Minimum number of MeSH terms")
    parser.add_argument("--sep-list", dest="separators", default="|;,/", help="Separators for list fields")
    parser.add_argument("--unknown-mode", dest="unknown_mode", action="store_true", help="Enable unknown label for ambiguous cases")
    parser.add_argument(
        "--prefer-pubmed-epsilon",
        dest="prefer_pubmed_epsilon",
        type=float,
        default=0.0,
        help="Optional epsilon added to score when PubMed votes review",
    )
    parser.add_argument("--chunksize", dest="chunksize", type=int, default=None, help="Optional chunk size for streaming input")
    parser.add_argument("--output-format", dest="output_format", choices=["csv", "parquet"], default="csv")
    parser.add_argument("--log-level", dest="log_level", default="INFO", help="Logging level")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level.upper(), logging.INFO))
    cfg_dict: Dict[str, Any] = {
        "separators": args.separators,
        "delta": args.delta,
        "k_min": args.k_min,
        "unknown_mode": args.unknown_mode,
        "prefer_pubmed_epsilon": args.prefer_pubmed_epsilon,
    }
    if args.config:
        with open(args.config, "r", encoding="utf-8") as fh:
            cfg_from_yaml = yaml.safe_load(fh) or {}
        cfg_dict.update(cfg_from_yaml)
    cfg = Config(**cfg_dict)
    LOGGER.info("Starting classification")
    run_pipeline(
        input_all=args.input_all,
        input_mesh=args.input_mesh,
        output=args.output,
        cfg=cfg,
        chunksize=args.chunksize,
        output_format=args.output_format,
    )
    LOGGER.info("Finished: output saved to %s", args.output)


if __name__ == "__main__":  # pragma: no cover
    main()
