"""Tests for the legacy artefact cleanup helpers."""

from __future__ import annotations

import json
from pathlib import Path

from library.maintenance import cleanup_legacy_outputs, ensure_legacy_cleanup


def _create_files(paths: list[Path]) -> None:
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("data", encoding="utf-8")


def test_cleanup_removes_known_patterns(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    keep_dataset = output_dir / "output.dataset.csv"
    keep_quality = output_dir / "output.dataset_quality_report_table.csv"
    legacy_json = output_dir / "output.dataset.quality.json"
    legacy_failures = output_dir / "output.dataset_failure_cases.csv"

    _create_files([keep_dataset, keep_quality, legacy_json, legacy_failures])

    result = cleanup_legacy_outputs(output_dir)

    assert result.removed_count == 2
    assert legacy_json.exists() is False
    assert legacy_failures.exists() is False
    assert keep_dataset.exists() is True
    assert keep_quality.exists() is True
    assert result.sentinel_path.exists() is True

    sentinel_payload = json.loads(result.sentinel_path.read_text(encoding="utf-8"))
    assert sentinel_payload["removed"]
    assert sentinel_payload["version"]


def test_cleanup_skips_when_sentinel_present(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    legacy_json = output_dir / "output.dataset.quality.json"
    _create_files([legacy_json])

    sentinel = output_dir / ".chembl-da-legacy-cleanup.json"
    sentinel.write_text("{}", encoding="utf-8")

    result = cleanup_legacy_outputs(output_dir)

    assert result.skipped is True
    assert legacy_json.exists() is True


def test_cleanup_dry_run_reports_without_deleting(tmp_path: Path) -> None:
    output_dir = tmp_path / "output"
    legacy_json = output_dir / "output.dataset.quality.json"
    _create_files([legacy_json])

    result = cleanup_legacy_outputs(output_dir, dry_run=True, force=True)

    assert result.dry_run is True
    assert result.removed_count == 1
    assert result.sentinel_path.exists() is False
    assert legacy_json.exists() is True


def test_ensure_legacy_cleanup_uses_config_output_dir(cfg) -> None:
    output_dir = Path(cfg.io.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    legacy_json = output_dir / "output.dataset.quality.json"
    legacy_json.write_text("legacy", encoding="utf-8")

    first_run = ensure_legacy_cleanup(cfg)
    assert first_run.removed_count == 1
    assert first_run.sentinel_path.exists() is True

    second_run = ensure_legacy_cleanup(cfg)
    assert second_run.skipped is True
