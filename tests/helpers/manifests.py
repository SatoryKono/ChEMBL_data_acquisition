"""Helpers for interacting with pipeline run manifests in tests."""

from __future__ import annotations

import json
from pathlib import Path


def list_manifest_files(base_path: Path) -> list[Path]:
    """Return all manifest files generated for the provided ``base_path``."""

    reports_dir = Path(base_path) / "reports"
    manifests = [
        path
        for path in reports_dir.glob("run_*.json")
        if path.name != "run_manifest.json"
    ]
    return sorted(manifests)


def load_latest_manifest(base_path: Path) -> tuple[Path, dict[str, object]]:
    """Load the most recent manifest ensuring the ``latest`` alias is consistent."""

    manifests = list_manifest_files(base_path)
    assert manifests, "expected run manifest to be written"
    latest = manifests[-1]

    alias = Path(base_path) / "reports" / "run_manifest.json"
    assert alias.exists(), "expected manifest alias to exist"
    if alias.is_symlink():
        assert (
            alias.resolve() == latest.resolve()
        ), "manifest alias must reference the latest run"
    else:
        assert alias.read_text(encoding="utf-8") == latest.read_text(
            encoding="utf-8"
        ), "manifest alias content must match the latest run"

    payload = json.loads(latest.read_text(encoding="utf-8"))
    return latest, payload
