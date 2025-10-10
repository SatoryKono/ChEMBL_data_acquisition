from __future__ import annotations

import csv
import hashlib
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest
import yaml

from library.pipelines.activity import get_activities
from scripts import get_activities as scripts_get_activities


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
def test_get_activities_cli_wrapper__script_entrypoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    limit = 4
    output_csv = tmp_path / "activities.csv"
    recorded: dict[str, object] = {}

    def fake_apply_config_overrides(
        args,
        parser,
        config_path,
        mapping=None,
        *,
        base_parser=None,
    ):
        del parser, config_path, mapping, base_parser
        args._config_metadata = None
        target = Path(args.final_out)
        args.final_out = target
        args.output_csv = target
        cfg = SimpleNamespace(
            io=SimpleNamespace(output_dir=target.parent),
            activity=SimpleNamespace(limit=None),
        )
        return cfg

    def fake_ensure_dirs(_cfg) -> None:
        return None

    def fake_run(_cfg, args) -> int:
        destination = Path(args.final_out)
        destination.parent.mkdir(parents=True, exist_ok=True)
        dataset = pd.DataFrame(
            [{"activity_id": f"ACT{i:03d}"} for i in range(limit)]
        )
        dataset.to_csv(destination, index=False, encoding="utf-8")
        quality_df = pd.DataFrame(
            [{"metric": "rows_total", "value": int(dataset.shape[0])}]
        )
        correlation_df = pd.DataFrame(
            [
                {
                    "column_a": "activity_id",
                    "column_b": "activity_id",
                    "correlation": 1.0,
                }
            ]
        )
        quality_path = destination.with_name(
            f"{destination.stem}_quality_report_table.csv"
        )
        correlation_path = destination.with_name(
            f"{destination.stem}_data_correlation_report_table.csv"
        )
        quality_df.to_csv(quality_path, index=False, encoding="utf-8")
        correlation_df.to_csv(
            correlation_path, index=False, encoding="utf-8"
        )
        recorded.update(
            {
                "dataset": dataset,
                "quality": quality_df,
                "correlation": correlation_df,
                "paths": {
                    "dataset": destination,
                    "quality": quality_path,
                    "correlation": correlation_path,
                },
            }
        )
        return 0

    monkeypatch.setattr(
        "library.utils.cli_tools.get_activities.cli.apply_config_overrides",
        fake_apply_config_overrides,
    )
    monkeypatch.setattr(
        "library.utils.cli_tools.get_activities.ensure_dirs", fake_ensure_dirs
    )
    monkeypatch.setattr(
        "library.utils.cli_tools.get_activities.run", fake_run
    )
    monkeypatch.setattr(
        "library.utils.cli_tools.get_activities.cli.configure_logger",
        lambda cfg: None,
    )

    exit_code = scripts_get_activities.main(
        [
            "--limit",
            str(limit),
            "--final-out",
            str(output_csv),
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()
    # Combine checks from both branches: validate output, metadata, and additional csv files.
    # 1. Validate main dataset CSV header
    with output_csv.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
    assert header.replace("\ufeff", "") == "activity_id"

    # Validate output file is not empty
    assert output_csv.stat().st_size > 0

    # 2. Validate metadata file existence and contents
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

    # 3. Validate that quality and correlation report files exist and are tracked
    assert recorded
    recorded_paths = recorded["paths"]
    assert recorded_paths["dataset"] == output_csv
    quality_path = recorded_paths["quality"]
    correlation_path = recorded_paths["correlation"]
    assert quality_path.exists()
    assert correlation_path.exists()
    assert {
        path.name for path in output_csv.parent.glob("*.csv")
    } == {
        output_csv.name,
        quality_path.name,
        correlation_path.name,
    }
