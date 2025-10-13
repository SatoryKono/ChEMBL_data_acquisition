"""Offline stub for the activity pipeline CLI used in integration tests."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def _read_fixture_rows(fixture_path: Path) -> list[list[str]]:
    with fixture_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        return list(reader)


def _write_limited_rows(
    destination: Path, rows: list[list[str]], limit: int | None
) -> None:
    header, *data_rows = rows
    if limit is not None and limit >= 0:
        data_rows = data_rows[:limit]
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(data_rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Lightweight offline implementation of the activity CLI."
    )
    parser.add_argument("--fixtures-dir", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=False)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args(argv)

    fixture_path = args.fixtures_dir / "activity.csv"
    if not fixture_path.exists():
        parser.error(f"Fixture file '{fixture_path}' is missing")

    rows = _read_fixture_rows(fixture_path)
    _write_limited_rows(args.destination, rows, args.limit)
    return 0


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    raise SystemExit(main())
