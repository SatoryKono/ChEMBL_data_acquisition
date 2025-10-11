"""CLI helper to purge legacy pipeline artefacts once after migration."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # pragma: no cover - import path guard
    sys.path.insert(0, str(ROOT))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Remove legacy metadata sidecars (`*.meta.yaml`, `.quality.json`, "
            "`*_failure_cases.csv`) created by historical runs. The cleanup "
            "runs automatically on first use after upgrading, but this helper "
            "allows manual execution or inspection."
        ),
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help=(
            "Configuration file used to resolve `io.output_dir`. When omitted "
            "the default search path is used."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override the output directory instead of loading it from config.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List files slated for removal without deleting or writing the sentinel.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Ignore the sentinel marker and re-run the cleanup.",
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="List matching legacy artefacts and exit without deleting them.",
    )
    return parser


def _resolve_output_dir(args: argparse.Namespace) -> Path:
    if isinstance(args.output_dir, Path):
        return args.output_dir

    config_path = args.config
    if config_path is None:
        from library.config import ConfigError, load_config

        try:
            cfg = load_config()
        except ConfigError as exc:
            print(f"Configuration error: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc
    else:
        from library.config import ConfigError, load_config

        try:
            cfg = load_config(config_path)
        except (ConfigError, FileNotFoundError) as exc:
            print(f"Configuration error for {config_path}: {exc}", file=sys.stderr)
            raise SystemExit(1) from exc

    return Path(cfg.io.output_dir)


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    from library.maintenance import cleanup_legacy_outputs, list_legacy_artifacts

    output_dir = _resolve_output_dir(args)

    if args.list_only:
        artefacts = list_legacy_artifacts(output_dir)
        for path in artefacts:
            print(path)
        print(f"{len(artefacts)} artefact(s) matched in {output_dir}")
        return 0

    result = cleanup_legacy_outputs(
        output_dir,
        dry_run=bool(args.dry_run),
        force=bool(args.force),
    )

    if result.skipped:
        reason = "sentinel present" if (result.sentinel_path.exists()) else "output directory missing"
        print(f"Cleanup skipped: {reason} ({result.output_dir})")
        return 0

    if result.dry_run:
        for path in result.removed:
            print(path)
        print(
            f"Dry run: {result.removed_count} artefact(s) would be removed from {result.output_dir}"
        )
        return 0

    if result.error_count:
        for path, message in result.errors:
            print(f"Failed to remove {path}: {message}", file=sys.stderr)
        print(
            f"Removed {result.removed_count} artefact(s); {result.error_count} error(s) encountered."
        )
        return 1

    print(
        "Removed {removed} legacy artefact(s) from {directory}. "
        "Future runs keep only the dataset and QA CSVs; re-enable diagnostics "
        "with --emit-legacy-artifacts/--debug/--keep-intermediate when needed."
        .format(removed=result.removed_count, directory=result.output_dir)
    )
    return 0


if __name__ == "__main__":  # pragma: no cover - manual execution entry point
    raise SystemExit(main())
