from __future__ import annotations

import argparse
import io
import json
import shutil
from collections import deque
from collections.abc import Callable
from pathlib import Path

import pandas as pd
import pytest

from library.common.csv_utils import sha256_file
from scripts import get_data


PipelineFunc = Callable[[list[str]], int]


def _build_stub_pipeline(
    name: str,
    key_column: str,
    required_columns: list[str],
    transform: Callable[[pd.DataFrame, get_data.Logger], pd.DataFrame],
    *,
    accept_subcommand: bool = False,
) -> PipelineFunc:
    """Return a stub ``main`` function emulating a pipeline step."""

    def _main(argv: list[str]) -> int:
        parser = argparse.ArgumentParser(prog=f"get_data_{name}_stub")
        if accept_subcommand:
            parser.add_argument("subcommand", nargs="?", default=None)
        parser.add_argument("--config", required=True)
        parser.add_argument("--input", required=True)
        parser.add_argument("--final-out", required=True)
        parser.add_argument("--log-level", default="INFO")
        parser.add_argument("--limit", type=int, default=None)
        parser.add_argument("--force", action="store_true")
        parser.add_argument("--skip-existing", action="store_true")
        args = parser.parse_args(argv)

        frame = pd.read_csv(Path(args.input))
        missing = [column for column in required_columns if column not in frame.columns]
        if missing:
            get_data._LOGGER.error(f"{name}_schema_invalid", missing=missing)
            return 1

        frame = frame[required_columns].copy()
        if key_column not in frame.columns:
            get_data._LOGGER.error(
                f"{name}_missing_key_column", column=key_column
            )
            return 1

        frame[key_column] = frame[key_column].astype("string").str.strip()
        duplicate_mask = frame[key_column].duplicated()
        if duplicate_mask.any():
            get_data._LOGGER.warning(
                f"{name}_duplicates_dropped", count=int(duplicate_mask.sum())
            )
            frame = frame.loc[~duplicate_mask].copy()

        transformed = transform(frame, get_data._LOGGER)
        transformed = transformed.sort_values(key_column).reset_index(drop=True)

        if transformed[key_column].duplicated().any():
            get_data._LOGGER.error(
                f"{name}_output_duplicates", count=int(transformed[key_column].duplicated().sum())
            )
            return 1

        Path(args.final_out).parent.mkdir(parents=True, exist_ok=True)
        transformed.to_csv(Path(args.final_out), index=False)
        return 0

    return _main


def _documents_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    pmids = pd.to_numeric(frame["pubmed_id"], errors="coerce").astype("Int64")
    missing_pubmed = pmids.isna()
    if int(missing_pubmed.sum()):
        logger.warning("document_missing_pubmed", count=int(missing_pubmed.sum()))
    frame = frame.copy()
    frame["title"] = frame["title"].astype("string").str.strip()
    frame["pubmed_id"] = pmids.astype("string").fillna("")
    frame["has_pubmed"] = (~missing_pubmed).astype("boolean")
    return frame[["document_chembl_id", "title", "pubmed_id", "has_pubmed"]]


def _targets_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    _ = logger  # logger kept for parity with other transforms
    frame = frame.copy()
    frame["target_name"] = frame["target_name"].astype("string").str.strip()
    frame["organism"] = frame["organism"].astype("string").str.strip()
    frame["name_upper"] = frame["target_name"].str.upper()
    return frame[["target_chembl_id", "target_name", "organism", "name_upper"]]


def _assays_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    _ = logger
    frame = frame.copy()
    frame["description"] = frame["description"].astype("string").str.strip()
    frame["description_length"] = frame["description"].str.len().astype("Int64")
    return frame[
        [
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
            "description_length",
        ]
    ]


def _testitems_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    names = frame["preferred_name"].astype("string").fillna("").str.strip()
    missing = names == ""
    if int(missing.sum()):
        logger.warning("testitem_missing_name", count=int(missing.sum()))
    normalized = names.str.lower().astype("string")
    normalized[missing] = pd.NA
    frame = frame.copy()
    frame["preferred_name"] = names.where(~missing, "")
    frame["normalized_name"] = normalized
    frame["is_named"] = (~missing).astype("boolean")
    return frame[["molecule_chembl_id", "preferred_name", "normalized_name", "is_named"]]


def _activities_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    values = frame["standard_value"].astype("string").str.strip()
    missing = values.isna() | (values == "")
    if int(missing.sum()):
        logger.error("activity_missing_value", count=int(missing.sum()))
    cleaned = values.fillna("")
    numeric = pd.to_numeric(cleaned.replace("", pd.NA), errors="coerce").astype("Float64")
    frame = frame.copy()
    frame["standard_value"] = cleaned
    frame["standard_units"] = frame["standard_units"].astype("string").fillna("").str.strip()
    frame["standard_value_numeric"] = numeric
    frame["is_valid"] = numeric.notna().astype("boolean")
    return frame[
        [
            "activity_id",
            "assay_chembl_id",
            "molecule_chembl_id",
            "standard_value",
            "standard_units",
            "standard_value_numeric",
            "is_valid",
        ]
    ]


@pytest.mark.e2e
def test_get_data_end_to_end__miniature_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    fixture_dir = Path(__file__).resolve().parents[1] / "data"
    for stem in ("document", "target", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{stem}.csv", input_dir / f"{stem}.csv")

    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n")

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))
    monkeypatch.setattr(get_data, "_warm_parent_catalog", lambda _cfg: None)

    log_streams: deque[io.StringIO] = deque()

    def _configure_logger_stub(cfg: get_data.LoggerConfig) -> get_data.Logger:
        stream = io.StringIO()
        log_streams.append(stream)
        return get_data.Logger(
            get_data.LoggerConfig(
                level=cfg.level,
                run_id=cfg.run_id,
                redact_secrets=cfg.redact_secrets,
                stream=stream,
            )
        )

    monkeypatch.setattr(get_data, "configure_logger", _configure_logger_stub)

    stub_steps = (
        get_data.PipelineStep(
            "document",
            _build_stub_pipeline(
                "document",
                "document_chembl_id",
                ["document_chembl_id", "title", "pubmed_id"],
                _documents_transform,
                accept_subcommand=True,
            ),
            "all",
        ),
        get_data.PipelineStep(
            "target",
            _build_stub_pipeline(
                "target",
                "target_chembl_id",
                ["target_chembl_id", "target_name", "organism"],
                _targets_transform,
                accept_subcommand=True,
            ),
            "all",
        ),
        get_data.PipelineStep(
            "assay",
            _build_stub_pipeline(
                "assay",
                "assay_chembl_id",
                [
                    "assay_chembl_id",
                    "target_chembl_id",
                    "document_chembl_id",
                    "description",
                ],
                _assays_transform,
            ),
            None,
        ),
        get_data.PipelineStep(
            "testitem",
            _build_stub_pipeline(
                "testitem",
                "molecule_chembl_id",
                ["molecule_chembl_id", "preferred_name"],
                _testitems_transform,
            ),
            None,
        ),
        get_data.PipelineStep(
            "activity",
            _build_stub_pipeline(
                "activity",
                "activity_id",
                [
                    "activity_id",
                    "assay_chembl_id",
                    "molecule_chembl_id",
                    "standard_value",
                    "standard_units",
                ],
                _activities_transform,
            ),
            None,
        ),
    )
    monkeypatch.setattr(get_data, "_PIPELINE_STEPS", stub_steps, raising=False)

    date_prefix = "20240102"
    def _invoke(argv: list[str]) -> tuple[int, list[dict[str, object]]]:
        exit_code = get_data.main(argv)
        assert log_streams, "expected logger stream to be captured"
        stream = log_streams.popleft()
        records = [
            json.loads(line)
            for line in stream.getvalue().splitlines()
            if line.strip()
        ]
        return exit_code, records

    argv = [
        "--base-path",
        str(base_path),
        "--input-dir",
        "input",
        "--output-dir",
        "output",
        "--config",
        str(config_path),
        "--date",
        date_prefix,
        "--log-level",
        "INFO",
        "--force",
    ]

    exit_code, logs = _invoke(argv)
    assert exit_code == 0

    events = {record.get("event"): record for record in logs if "event" in record}
    assert "document_duplicates_dropped" in events
    assert events["document_duplicates_dropped"].get("level") == "WARN"
    assert "activity_missing_value" in events
    assert events["activity_missing_value"].get("level") == "ERROR"

    expected_dir = Path(__file__).resolve().parents[1] / "resources" / "expected_get_data"
    key_columns = {
        "document": "document_chembl_id",
        "target": "target_chembl_id",
        "assay": "assay_chembl_id",
        "testitem": "molecule_chembl_id",
        "activity": "activity_id",
    }
    output_paths: dict[str, Path] = {}
    for step_name, stem in get_data._DEFAULT_OUTPUT_STEMS.items():
        final_path = output_dir / f"output.{stem}_{date_prefix}.csv"
        output_paths[step_name] = final_path
        assert final_path.exists(), f"expected output for {step_name} missing"
        actual = pd.read_csv(final_path)
        expected = pd.read_csv(expected_dir / f"{stem}.csv")
        pd.testing.assert_frame_equal(actual, expected)
        assert not actual[key_columns[step_name]].duplicated().any()

    hashes_before = {name: sha256_file(path) for name, path in output_paths.items()}

    repeat_exit_code, repeat_logs = _invoke(argv)
    assert repeat_exit_code == 0
    for step_name in key_columns:
        assert any(
            record["event"] == "step_done" and record.get("step") == step_name
            for record in repeat_logs
        )

    hashes_after = {name: sha256_file(path) for name, path in output_paths.items()}
    assert hashes_before == hashes_after

    degraded_input_dir = base_path / "degraded_input"
    degraded_output_dir = base_path / "degraded_output"
    shutil.copytree(input_dir, degraded_input_dir)
    degraded_output_dir.mkdir()

    degraded_document = degraded_input_dir / "document.csv"
    broken_document = pd.read_csv(degraded_document).drop(columns=["pubmed_id"])
    broken_document.to_csv(degraded_document, index=False)

    degraded_argv = [
        "--base-path",
        str(base_path),
        "--input-dir",
        "degraded_input",
        "--output-dir",
        "degraded_output",
        "--config",
        str(config_path),
        "--date",
        date_prefix,
        "--log-level",
        "INFO",
        "--force",
    ]

    degraded_exit_code, degraded_logs = _invoke(degraded_argv)
    assert degraded_exit_code == 1
    degraded_events = {record.get("event"): record for record in degraded_logs if "event" in record}
    assert "document_schema_invalid" in degraded_events
    assert degraded_events["document_schema_invalid"].get("level") == "ERROR"
    assert "workflow_failed" in degraded_events
    assert not any(degraded_output_dir.glob("*.csv"))
