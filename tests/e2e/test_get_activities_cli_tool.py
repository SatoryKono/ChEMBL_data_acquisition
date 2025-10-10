from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from library.pipelines.activity import get_activities


@pytest.mark.e2e
def test_get_activities_cli__subprocess_execution(
    tmp_path: Path, record_property
) -> None:
    limit = 5
    project_root = Path(__file__).resolve().parents[2]

    env = os.environ.copy()
    pythonpath_parts = [str(project_root)]
    existing = env.get("PYTHONPATH")
    if existing:
        pythonpath_parts.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    command = [
        sys.executable,
        "-m",
        "library.utils.cli_tools.get_activities",
        "--limit",
        str(limit),
    ]

    completed = subprocess.run(
        command,
        cwd=project_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert (
        completed.returncode == 0
    ), f"CLI execution failed: {completed.stderr or completed.stdout}"

    stdout_lines = [line for line in completed.stdout.splitlines() if line.strip()]
    assert any(
        "generated" in line and f"count={limit}" in line for line in stdout_lines
    ), completed.stdout

    rows = get_activities(limit)
    csv_path = tmp_path / "activities.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["activity_id"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    expected_hash = "8ab1f75124c8cac1102f311897e36d9902a947b999b7143536281c266d98b28e"
    assert digest == expected_hash

    record_property("artifact_csv", str(csv_path))
    record_property("artifact_hash", digest)


@pytest.mark.e2e
def test_get_activities_cli_wrapper__script_entrypoint(tmp_path: Path) -> None:
    limit = 4
    project_root = Path(__file__).resolve().parents[2]

    env = os.environ.copy()
    pythonpath_parts = [str(project_root)]
    existing = env.get("PYTHONPATH")
    if existing:
        pythonpath_parts.append(existing)
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_parts)

    output_csv = tmp_path / "activities.csv"
    command = [
        sys.executable,
        "scripts/get_activities.py",
        "--limit",
        str(limit),
        "--final-out",
        str(output_csv),
    ]

    completed = subprocess.run(
        command,
        cwd=project_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    assert (
        completed.returncode == 0
    ), f"CLI wrapper execution failed: {completed.stderr or completed.stdout}"

    assert output_csv.exists()
    with output_csv.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
    assert header.replace("\ufeff", "") == "activity_id"

    meta_path = output_csv.with_suffix(output_csv.suffix + ".meta.yaml")
    assert meta_path.exists()
    metadata = yaml.safe_load(meta_path.read_text(encoding="utf-8"))
    assert isinstance(metadata, dict)
    base_keys = {
        "generated_at",
        "git_sha",
        "python_version",
        "platform",
        "command",
        "config",
        "inputs",
        "stats",
        "schema",
        "columns",
        "dtypes",
        "pipeline_version",
    }
    optional_keys = {"invocation", "dictionaries"}
    assert base_keys.issubset(metadata)
    assert set(metadata).issubset(base_keys | optional_keys)
    assert metadata.get("schema") in {None, "ActivitiesSchema"}
    assert isinstance(metadata.get("stats"), dict)
