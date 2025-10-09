#!/usr/bin/env python3
"""Verify deterministic CSV output.

This script writes a small test :class:`pandas.DataFrame` using both the
standard and chunked deterministic CSV writers and compares the SHA-256 hashes
of the resulting files. A non-zero exit code is returned if any hash differs.
"""

from __future__ import annotations

# ruff: noqa: E402
import argparse
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

try:
    import pandas as pd
except ImportError as exc:  # pragma: no cover - import-time check
    raise SystemExit(
        "pandas is required to run library.utils.cli_tools.check_determinism."
        " Install it with 'pip install pandas'."
    ) from exc

from library.cli import configure_logger, create_logger_config
from library.cli.base import compute_generated_at
from library.common.csv_utils import (
    sha256_file,
    write_csv_chunks_deterministic,
    write_csv_deterministic,
)
from library.common.log import logger
from library.common.timing import log_duration


def run_check(tmp_dir: Path) -> bool:
    """Write the same CSV twice and compare their SHA-256 hashes.

    Parameters
    ----------
    tmp_dir:
        Directory where temporary CSV files will be created.

    Returns
    -------
    bool
        ``True`` if both files produce identical hashes, ``False`` otherwise.
    """

    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})

    first = tmp_dir / "first.csv"
    second = tmp_dir / "second.csv"
    chunked = tmp_dir / "chunked.csv"

    # Write the DataFrame twice using the deterministic writer
    key_cols = list(df.columns)
    write_csv_deterministic(df, first, key_cols=key_cols)
    hash1 = sha256_file(first)

    write_csv_deterministic(df, second, key_cols=key_cols)
    hash2 = sha256_file(second)

    chunk_iter = [df.copy()]
    write_csv_chunks_deterministic(
        chunk_iter,
        chunked,
        key_cols=key_cols,
        col_order=key_cols,
        chunksize=len(df),
        merge_chunksize=len(df),
        sort_chunksize=len(df),
        sep=",",
        encoding="utf-8-sig",
        cfg=None,
    )
    hash3 = sha256_file(chunked)

    logger.debug("hash", label="first", value=hash1)
    logger.debug("hash", label="second", value=hash2)
    logger.debug("hash", label="chunked", value=hash3)

    metadata_status = "skipped"
    first_meta = first.with_suffix(first.suffix + ".meta.yaml")
    second_meta = second.with_suffix(second.suffix + ".meta.yaml")
    chunked_meta = chunked.with_suffix(chunked.suffix + ".meta.yaml")

    meta_entries: list[tuple[str, Path, str]] = []
    for label, path in (
        ("first_meta", first_meta),
        ("second_meta", second_meta),
        ("chunked_meta", chunked_meta),
    ):
        if path.exists():
            digest = sha256_file(path)
            logger.debug("hash", label=label, value=digest)
            meta_entries.append((label, path, digest))

    if len(meta_entries) >= 2:
        metadata_status = "matched"
        baseline = meta_entries[0][2]
        if any(digest != baseline for _, _, digest in meta_entries[1:]):
            metadata_status = "mismatch"
            logger.warning(
                "metadata_hash_mismatch",
                entries=[
                    {"label": label, "path": str(path), "hash": digest}
                    for label, path, digest in meta_entries
                ],
            )
            logger.info("metadata_check", status=metadata_status)
            return False

    logger.info("metadata_check", status=metadata_status)

    return hash1 == hash2 == hash3


def main() -> int:
    """CLI entry point.

    Returns
    -------
    int
        ``0`` on success, ``1`` when hashes differ.
    """

    start = perf_counter()
    try:
        parser = argparse.ArgumentParser(
            description="Check deterministic CSV writing",
        )
        parser.add_argument(
            "--log-level",
            default="INFO",
            help="Logging level (default: INFO).",
        )
        args = parser.parse_args()

        log_cfg = create_logger_config(args.log_level)
        seed_parts = [parser.prog]
        seed_parts.extend(sys.argv[1:])
        log_cfg.generated_at = compute_generated_at(
            date_token=None,
            run_id=log_cfg.run_id,
            seed_parts=seed_parts,
        )
        configure_logger(log_cfg)

        with TemporaryDirectory() as tmp:
            ok = run_check(Path(tmp))

        if ok:
            logger.info("hashes_match")
            return 0

        logger.error("hashes_differ")
        return 1
    finally:
        log_duration(start)


if __name__ == "__main__":
    raise SystemExit(main())
