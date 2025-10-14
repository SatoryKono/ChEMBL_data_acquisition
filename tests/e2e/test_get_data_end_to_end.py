from __future__ import annotations

import argparse
import io
import json
import re
import shutil
import textwrap
from collections import deque
from collections.abc import Callable
from functools import lru_cache
from pathlib import Path

import pandas as pd
import pytest

from library.common.csv_utils import sha256_file
from library.project_version import get_pipeline_version
from scripts import get_data
from tests.helpers import ASSAY_ENRICHMENT_MIN_RATIO
from tests.helpers.logs import parse_log_file, parse_log_lines

PipelineFunc = Callable[[list[str]], int]


E2E_TEST_GIT_SHA = "e2e-test-sha"


def _build_stub_pipeline(
    name: str,
    key_column: str,
    required_columns: list[str],
    transform: Callable[[pd.DataFrame, get_data.Logger], pd.DataFrame],
    *,
    accept_subcommand: bool = False,
    accept_mode: bool = False,
    optional_columns: list[str] | None = None,
    post_process: Callable[[Path], None] | None = None,
) -> PipelineFunc:
    """Return a stub ``main`` function emulating a pipeline step."""

    def _main(argv: list[str]) -> int:
        parser = argparse.ArgumentParser(prog=f"get_data_{name}_stub")
        if accept_subcommand:
            parser.add_argument("subcommand", nargs="?", default=None)
        if accept_mode:
            parser.add_argument("--mode", required=True)
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

        selected_columns = list(required_columns)
        if optional_columns:
            selected_columns.extend(
                column for column in optional_columns if column in frame.columns
            )
        frame = frame[selected_columns].copy()
        if key_column not in frame.columns:
            get_data._LOGGER.error(f"{name}_missing_key_column", column=key_column)
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
                f"{name}_output_duplicates",
                count=int(transformed[key_column].duplicated().sum()),
            )
            return 1

        destination = Path(args.final_out)
        destination.parent.mkdir(parents=True, exist_ok=True)
        transformed.to_csv(destination, index=False)
        if post_process is not None:
            post_process(destination)
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
    lookup = _load_assay_dictionary()
    enriched = frame.merge(lookup, on="assay_chembl_id", how="left")
    quality_columns = ["assay_strain", "assay_group", "year", "accession"]
    completeness = 1.0 - enriched[quality_columns].isna().mean()
    min_ratio = float(completeness.min()) if len(completeness) else 1.0
    if min_ratio < ASSAY_ENRICHMENT_MIN_RATIO:
        raise AssertionError(
            "assay enrichment below threshold "
            f"(threshold={ASSAY_ENRICHMENT_MIN_RATIO}, completeness={completeness.to_dict()})"
        )
    enriched["year"] = enriched["year"].astype("Int64")
    return enriched[
        [
            "assay_chembl_id",
            "target_chembl_id",
            "document_chembl_id",
            "description",
            "description_length",
            "assay_strain",
            "assay_group",
            "year",
            "accession",
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
    return frame[
        ["molecule_chembl_id", "preferred_name", "normalized_name", "is_named"]
    ]


def _activities_transform(frame: pd.DataFrame, logger: get_data.Logger) -> pd.DataFrame:
    values = frame["standard_value"].astype("string").str.strip()
    missing = values.isna() | (values == "")
    if int(missing.sum()):
        logger.error("activity_missing_value", count=int(missing.sum()))
    if "force_failure" in frame.columns and frame["force_failure"].any():
        raise RuntimeError("activity_forced_failure")
    cleaned = values.fillna("")
    numeric = pd.to_numeric(cleaned.replace("", pd.NA), errors="coerce").astype(
        "Float64"
    )
    frame = frame.copy()
    frame["standard_value"] = cleaned
    frame["standard_units"] = (
        frame["standard_units"].astype("string").fillna("").str.strip()
    )
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


def _failing_target_transform(
    frame: pd.DataFrame, logger: get_data.Logger
) -> pd.DataFrame:
    _ = frame, logger
    raise RuntimeError("target_forced_failure")


@pytest.mark.e2e
@pytest.mark.smoke
def test_get_data_end_to_end__miniature_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    request: pytest.FixtureRequest,
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    cleanup_artifact = output_dir / "sentinel__cleanup.tmp"
    cleanup_artifact.write_text("cleanup", encoding="utf-8")
    external_dir = base_path / "external"
    external_dir.mkdir()
    external_artifact = external_dir / "sentinel__external.tmp"
    external_artifact.write_text("outside", encoding="utf-8")

    fixture_dir = Path(__file__).resolve().parents[1] / "resources" / "pipeline_inputs"
    for stem in ("document", "target", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{stem}.csv", input_dir / f"{stem}.csv")

    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n")

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))
    monkeypatch.setenv("GIT_SHA", E2E_TEST_GIT_SHA)
    monkeypatch.setattr(get_data, "_warm_parent_catalog", lambda _cfg, _base: None)

    log_streams: deque[io.StringIO] = deque()
    observed_run_ids: list[str] = []

    original_configure = get_data.configure_logger

    def _configure_logger_stub(cfg: get_data.LoggerConfig) -> get_data.Logger:
        stream = io.StringIO()
        log_streams.append(stream)
        observed_run_ids.append(str(cfg.run_id))
        return original_configure(
            get_data.LoggerConfig(
                level=cfg.level,
                run_id=cfg.run_id,
                redact_secrets=cfg.redact_secrets,
                stream=stream,
                handlers=list(cfg.handlers) if cfg.handlers else [],
                logger_name=cfg.logger_name,
            )
        )

    monkeypatch.setattr(get_data, "configure_logger", _configure_logger_stub)

    def _make_report_writer(step_count: int) -> Callable[[Path], None]:
        def _writer(working_output: Path) -> None:
            base_dir = working_output.parent.parent
            reports_dir = base_dir / "reports"
            reports_dir.mkdir(parents=True, exist_ok=True)
            run_id = observed_run_ids[-1] if observed_run_ids else "unknown"
            timestamp = "2020-01-01T00:00:00+00:00"
            payload = {
                "meta": {
                    "repo": "SatoryKono/ChEMBL_data_acquisition",
                    "commit": "0000000",
                    "branch": "test",
                    "ts_utc": timestamp,
                    "run_id": run_id,
                    "python": "3.11",
                    "pytest": "7.0",
                    "duration_sec": 0.0,
                },
                "summary": {
                    "total": step_count,
                    "passed": step_count,
                    "failed": 0,
                    "skipped": 0,
                    "xfailed": 0,
                    "xpassed": 0,
                    "error": 0,
                    "success_rate": 1.0,
                },
                "tests": [],
            }
            (reports_dir / "test_report.json").write_text(
                json.dumps(payload, indent=2), encoding="utf-8"
            )
            summary_md = (
                "# Test Summary\n\n"
                "- Repo: `SatoryKono/ChEMBL_data_acquisition`\n"
                "- Commit: 0000000\n"
                "- Branch: test\n"
                f"- Timestamp (UTC): {timestamp}\n"
                "- Duration: 0.0 s\n"
                "- Success rate: 100.00%\n\n"
                "| total | passed | failed | skipped | xfailed | xpassed | error |\n"
                "|------:|-------:|-------:|--------:|--------:|--------:|------:|\n"
                f"|{step_count:6d}|{step_count:7d}|{0:7d}|{0:8d}|{0:8d}|{0:8d}|{0:6d}|\n"
                "\n## Failed / Error details\n"
            )
            (reports_dir / "test_summary.md").write_text(summary_md, encoding="utf-8")

        return _writer

    report_writer = _make_report_writer(5)
    stub_steps = (
        get_data.PipelineStep(
            name="document",
            main=_build_stub_pipeline(
                "document",
                "document_chembl_id",
                ["document_chembl_id", "title", "pubmed_id"],
                _documents_transform,
                accept_mode=True,
            ),
            input_filename="document.csv",
            output_stem="documents",
            extra_args=("--mode", "all"),
        ),
        get_data.PipelineStep(
            name="target",
            main=_build_stub_pipeline(
                "target",
                "target_chembl_id",
                ["target_chembl_id", "target_name", "organism"],
                _targets_transform,
                accept_subcommand=True,
            ),
            input_filename="target.csv",
            output_stem="targets",
            subcommand="all",
        ),
        get_data.PipelineStep(
            name="assay",
            main=_build_stub_pipeline(
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
            input_filename="assay.csv",
            output_stem="assays",
        ),
        get_data.PipelineStep(
            name="testitem",
            main=_build_stub_pipeline(
                "testitem",
                "molecule_chembl_id",
                ["molecule_chembl_id", "preferred_name"],
                _testitems_transform,
            ),
            input_filename="testitem.csv",
            output_stem="testitem",
        ),
        get_data.PipelineStep(
            name="activity",
            main=_build_stub_pipeline(
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
                optional_columns=["force_failure"],
                post_process=report_writer,
            ),
            input_filename="activity.csv",
            output_stem="activities",
        ),
    )
    monkeypatch.setattr(
        get_data,
        "_resolve_pipeline_steps",
        lambda _: stub_steps,
    )
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {}, raising=False)

    date_prefix = "20240102"

    def _invoke(argv: list[str]) -> tuple[int, list[dict[str, object]], str]:
        before = len(observed_run_ids)
        exit_code = get_data.main(argv)
        assert log_streams, "expected logger stream to be captured"
        records: list[dict[str, object]] = []
        while log_streams:
            stream = log_streams.popleft()
            records.extend(parse_log_lines(stream.getvalue()))
        assert len(observed_run_ids) > before, "run identifier was not captured"
        run_id = observed_run_ids[-1]
        return exit_code, records, run_id

    def _collect_log_run_ids(entries: list[dict[str, object]]) -> set[str]:
        collected: set[str] = set()
        for entry in entries:
            data = entry.get("data", {})
            if not isinstance(data, dict):
                continue
            value = data.get("run_id")
            if value is None:
                continue
            collected.add(str(value))
        return collected

    def _read_manifest_run_id() -> str:
        manifest_alias = base_path / "reports" / "run_manifest.json"
        assert manifest_alias.exists(), "expected manifest alias to exist"
        payload = json.loads(manifest_alias.read_text(encoding="utf-8"))
        run_info = payload["run"]
        assert run_info["pipeline_version"] == get_pipeline_version()
        assert run_info["git_sha"] == E2E_TEST_GIT_SHA
        value = payload["run"].get("run_id")
        assert value, "manifest must include run_id"
        return str(value)

    sentinel_path = get_data.OUTPUT_DIR / "cleanup_sentinel.tmp"
    sentinel_path.parent.mkdir(parents=True, exist_ok=True)
    sentinel_path.write_text("keep", encoding="utf-8")
    request.addfinalizer(lambda: sentinel_path.unlink(missing_ok=True))

    cleanup_file = output_dir / "temporary.tmp"
    cleanup_file.write_text("stale", encoding="utf-8")
    request.addfinalizer(lambda: cleanup_file.unlink(missing_ok=True))

    cleanup_dir = output_dir / "raw"
    cleanup_dir.mkdir(exist_ok=True)
    (cleanup_dir / "scratch.csv").write_text("id\n1\n", encoding="utf-8")
    request.addfinalizer(lambda: shutil.rmtree(cleanup_dir, ignore_errors=True))

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

    exit_code, logs, first_run_id = _invoke(argv)
    assert exit_code == 0
    assert _collect_log_run_ids(logs) == {first_run_id}

    assert not cleanup_file.exists()
    assert not cleanup_dir.exists()
    assert sentinel_path.exists()

    manifest_run_id = _read_manifest_run_id()
    assert manifest_run_id == first_run_id

    log_dir = base_path / "logs"
    orchestrator_log = log_dir / f"get_data_{date_prefix}.log"
    assert orchestrator_log.exists()
    orchestrator_records = parse_log_file(orchestrator_log)
    orchestrator_events = {entry.get("event") for entry in orchestrator_records}
    assert "pipeline_start" in orchestrator_events
    assert "workflow_succeeded" in orchestrator_events

    csv_summary_entries = [
        entry
        for entry in orchestrator_records
        if "Найдено" in entry.get("raw", "")
    ]
    assert csv_summary_entries, "expected CSV summary entry in orchestrator log"
    summary_message = csv_summary_entries[-1]["raw"].split(" | ", 2)[-1]
    resolved_output = str(output_dir.resolve())
    assert resolved_output in summary_message
    count_match = re.search(r"Найдено\s+(?P<count>\d+)\s+CSV-файлов", summary_message)
    assert count_match, f"unable to parse CSV count from '{summary_message}'"
    reported_count = int(count_match.group("count"))
    actual_csv_count = sum(1 for path in output_dir.glob("*.csv"))
    assert reported_count == actual_csv_count

    events = {record.get("event"): record for record in logs if "event" in record}
    assert "document_duplicates_dropped" in events
    assert events["document_duplicates_dropped"].get("level") == "WARNING"
    assert "activity_missing_value" in events
    assert events["activity_missing_value"].get("level") == "ERROR"
    assert not any(record.get("event") == "step_arguments" for record in logs)

    expected_dir = (
        Path(__file__).resolve().parents[1] / "resources" / "expected_get_data"
    )
    key_columns = {
        "document": "document_chembl_id",
        "target": "target_chembl_id",
        "assay": "assay_chembl_id",
        "testitem": "molecule_chembl_id",
        "activity": "activity_id",
    }
    output_paths: dict[str, Path] = {}
    for step_name, stem in get_data.DEFAULT_OUTPUT_STEMS.items():
        final_path = output_dir / f"output.{stem}_{date_prefix}.csv"
        output_paths[step_name] = final_path
        assert final_path.exists(), f"expected output for {step_name} missing"
        actual = pd.read_csv(final_path)
        expected = pd.read_csv(expected_dir / f"{stem}.csv")
        pd.testing.assert_frame_equal(actual, expected)
        assert not actual[key_columns[step_name]].duplicated().any()
        if step_name == "assay":
            quality_columns = ["assay_strain", "assay_group", "year", "accession"]
            completeness = 1.0 - actual[quality_columns].isna().mean()
            assert (completeness >= ASSAY_ENRICHMENT_MIN_RATIO).all(), completeness

    actual_csv_files = sorted(output_dir.glob("*.csv"))
    assert reported_csv_count == len(actual_csv_files)

    hashes_before = {name: sha256_file(path) for name, path in output_paths.items()}

    repeat_exit_code, repeat_logs, repeat_run_id = _invoke(argv)
    assert repeat_exit_code == 0
    assert repeat_run_id == first_run_id
    assert _collect_log_run_ids(repeat_logs) == {repeat_run_id}
    assert _read_manifest_run_id() == repeat_run_id
    for step_name in key_columns:
        assert any(
            record["event"] == "step_done"
            and record.get("data", {}).get("step") == step_name
            for record in repeat_logs
        )
    assert not any(record.get("event") == "step_arguments" for record in repeat_logs)

    hashes_after = {name: sha256_file(path) for name, path in output_paths.items()}
    assert hashes_before == hashes_after

    verbose_argv = [*argv, "--verbose"]
    verbose_exit_code, verbose_logs, verbose_run_id = _invoke(verbose_argv)
    assert verbose_exit_code == 0
    assert any(record.get("event") == "step_arguments" for record in verbose_logs)
    assert any(record.get("level") == "DEBUG" for record in verbose_logs)
    assert _collect_log_run_ids(verbose_logs) == {verbose_run_id}
    verbose_hashes = {name: sha256_file(path) for name, path in output_paths.items()}
    assert hashes_before == verbose_hashes
    verbose_file_records = parse_log_file(orchestrator_log)
    assert any(entry.get("event") == "step_arguments" for entry in verbose_file_records)
    assert _collect_log_run_ids(verbose_file_records) == {first_run_id, verbose_run_id}
    assert _read_manifest_run_id() == verbose_run_id

    reports_dir = base_path / "reports"
    report_json_path = reports_dir / "test_report.json"
    summary_md_path = reports_dir / "test_summary.md"
    assert report_json_path.exists()
    assert summary_md_path.exists()
    report_payload = json.loads(report_json_path.read_text(encoding="utf-8"))
    assert report_payload["meta"]["repo"] == "SatoryKono/ChEMBL_data_acquisition"
    assert report_payload["meta"]["run_id"] == verbose_run_id
    assert report_payload["meta"]["ts_utc"] == "2020-01-01T00:00:00+00:00"
    assert report_payload["summary"]["total"] == 5
    assert report_payload["summary"]["passed"] == 5
    summary_md = summary_md_path.read_text(encoding="utf-8")
    assert "Success rate: 100.00%" in summary_md
    assert "|     5|      5|" in summary_md

    new_date_prefix = "20240103"
    new_run_argv = [
        "--base-path",
        str(base_path),
        "--input-dir",
        "input",
        "--output-dir",
        "output",
        "--config",
        str(config_path),
        "--date",
        new_date_prefix,
        "--log-level",
        "INFO",
        "--force",
    ]

    new_exit_code, _, new_run_id = _invoke(new_run_argv)
    assert new_exit_code == 0
    assert new_run_id != verbose_run_id
    new_log_path = log_dir / f"get_data_{new_date_prefix}.log"
    assert new_log_path.exists()
    for _step_name, stem in get_data.DEFAULT_OUTPUT_STEMS.items():
        new_path = output_dir / f"output.{stem}_{new_date_prefix}.csv"
        assert new_path.exists()

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

    degraded_exit_code, degraded_logs, degraded_run_id = _invoke(degraded_argv)
    assert degraded_exit_code == 1
    assert degraded_run_id not in {first_run_id, verbose_run_id, new_run_id}
    degraded_events = {
        record.get("event"): record for record in degraded_logs if "event" in record
    }
    assert "document_schema_invalid" in degraded_events
    assert degraded_events["document_schema_invalid"].get("level") == "ERROR"
    assert "workflow_failed" in degraded_events
    assert not any(degraded_output_dir.glob("*.csv"))

    partial_base = base_path / "partial"
    partial_input = partial_base / "input"
    partial_output = partial_base / "output"
    partial_input.mkdir(parents=True)
    partial_output.mkdir()
    for stem in ("document", "target", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{stem}.csv", partial_input / f"{stem}.csv")
    activity_partial = pd.read_csv(partial_input / "activity.csv")
    activity_partial["force_failure"] = True
    activity_partial.to_csv(partial_input / "activity.csv", index=False)
    partial_config = partial_base / "config.yaml"
    partial_config.write_text("io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n")
    partial_argv = [
        "--base-path",
        str(partial_base),
        "--input-dir",
        "input",
        "--output-dir",
        "output",
        "--config",
        str(partial_config),
        "--date",
        date_prefix,
        "--log-level",
        "INFO",
        "--force",
    ]
    partial_exit_code, partial_logs, partial_run_id = _invoke(partial_argv)
    assert partial_exit_code == 1
    assert partial_run_id not in {
        first_run_id,
        repeat_run_id,
        verbose_run_id,
        new_run_id,
        degraded_run_id,
    }
    partial_events = {
        record.get("event"): record for record in partial_logs if "event" in record
    }
    assert "step_exception" in partial_events
    assert (partial_output / f"output.documents_{date_prefix}.csv").exists()
    activity_sentinel = partial_output / f"output.activity_{date_prefix}.csv.failed"
    assert activity_sentinel.exists()

    missing_base = base_path / "missing"
    missing_input = missing_base / "input"
    missing_output = missing_base / "output"
    missing_input.mkdir(parents=True)
    missing_output.mkdir()
    for stem in ("document", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{stem}.csv", missing_input / f"{stem}.csv")
    config_missing = missing_base / "config.yaml"
    config_missing.write_text("io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n")
    missing_argv = [
        "--base-path",
        str(missing_base),
        "--input-dir",
        "input",
        "--output-dir",
        "output",
        "--config",
        str(config_missing),
        "--date",
        date_prefix,
        "--log-level",
        "INFO",
    ]
    missing_exit_code, missing_logs, missing_run_id = _invoke(missing_argv)
    assert missing_exit_code == 1
    assert missing_run_id not in {
        first_run_id,
        repeat_run_id,
        verbose_run_id,
        new_run_id,
        degraded_run_id,
        partial_run_id,
    }
    missing_events = {
        record.get("event"): record for record in missing_logs if "event" in record
    }
    assert missing_events.get("step_input_missing") is not None
    missing_target = missing_output / f"output.targets_{date_prefix}.csv"
    assert not missing_target.exists()
    sentinel_path = missing_output / f"output.targets_{date_prefix}.csv.failed"
    assert sentinel_path.exists()


@pytest.mark.e2e
def test_get_data_end_to_end__blocked_steps_after_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    base_path = tmp_path
    input_dir = base_path / "input"
    output_dir = base_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()

    fixture_dir = Path(__file__).resolve().parents[1] / "resources" / "pipeline_inputs"
    for stem in ("document", "target", "assay", "testitem", "activity"):
        shutil.copy(fixture_dir / f"{stem}.csv", input_dir / f"{stem}.csv")

    config_path = base_path / "config.yaml"
    config_path.write_text("io:\n  csv_sep: ','\n  csv_encoding: 'utf-8'\n")

    monkeypatch.setenv("CHEMBL_DA_BASE_PATH", str(base_path))
    monkeypatch.setattr(get_data, "_warm_parent_catalog", lambda _cfg, _base: None)
    monkeypatch.setenv("GIT_SHA", E2E_TEST_GIT_SHA)

    stub_steps = (
        get_data.PipelineStep(
            name="document",
            main=_build_stub_pipeline(
                "document",
                "document_chembl_id",
                ["document_chembl_id", "title", "pubmed_id"],
                _documents_transform,
                accept_mode=True,
            ),
            input_filename="document.csv",
            output_stem="documents",
            extra_args=("--mode", "all"),
        ),
        get_data.PipelineStep(
            name="target",
            main=_build_stub_pipeline(
                "target",
                "target_chembl_id",
                ["target_chembl_id", "target_name", "organism"],
                _failing_target_transform,
                accept_subcommand=True,
            ),
            input_filename="target.csv",
            output_stem="targets",
            subcommand="all",
        ),
        get_data.PipelineStep(
            name="assay",
            main=_build_stub_pipeline(
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
            input_filename="assay.csv",
            output_stem="assays",
        ),
        get_data.PipelineStep(
            name="testitem",
            main=_build_stub_pipeline(
                "testitem",
                "molecule_chembl_id",
                ["molecule_chembl_id", "preferred_name"],
                _testitems_transform,
            ),
            input_filename="testitem.csv",
            output_stem="testitem",
        ),
        get_data.PipelineStep(
            name="activity",
            main=_build_stub_pipeline(
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
            input_filename="activity.csv",
            output_stem="activities",
        ),
    )

    monkeypatch.setattr(
        get_data,
        "_resolve_pipeline_steps",
        lambda _: stub_steps,
    )
    monkeypatch.setattr(get_data, "_PIPELINE_APIS", {}, raising=False)

    date_prefix = "20240105"
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

    exit_code = get_data.main(argv)
    assert exit_code == 1

    manifest_path = base_path / "reports" / "run_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    run_info = manifest["run"]
    assert run_info["pipeline_version"] == get_pipeline_version()
    assert run_info["git_sha"] == E2E_TEST_GIT_SHA
    entries_by_name = {entry["name"]: entry for entry in manifest["steps"]}

    assert entries_by_name["document"]["status"] == "success"
    assert entries_by_name["target"]["status"] == "failed"
    assert entries_by_name["target"]["reason"] == "exception"

    for blocked_name in ("assay", "testitem", "activity"):
        blocked_entry = entries_by_name[blocked_name]
        assert blocked_entry["status"] == "blocked"
        assert blocked_entry["executed"] is False
        assert blocked_entry["reason"] == "dependency_failed"
        assert blocked_entry["exit_code"] is None

        final_output = (
            output_dir
            / f"output.{get_data.DEFAULT_OUTPUT_STEMS[blocked_name]}_{date_prefix}.csv"
        )
        assert not final_output.exists()

    target_sentinel = (
        output_dir
        / f"output.{get_data.DEFAULT_OUTPUT_STEMS['target']}_{date_prefix}.csv.failed"
    )
    assert target_sentinel.exists()


@pytest.mark.e2e
@pytest.mark.pipeline_scenario("enrichment")
def test_get_data_end_to_end__pubchem_disable_flag(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    stub_script = tmp_path / "stub_testitem_script.py"
    stub_script.write_text(
        textwrap.dedent(
            """
            from __future__ import annotations

            import argparse
            import json
            import os
            import sys
            from pathlib import Path


            def _pc_init_session(*_args: object, **_kwargs: object) -> None:
                raise SystemExit("pc.init_session should not be called when disabled")


            def main(argv: list[str] | None = None) -> int:
                parser = argparse.ArgumentParser()
                parser.add_argument(
                    "--pubchem-enable",
                    dest="pubchem_enable",
                    action=argparse.BooleanOptionalAction,
                    default=None,
                )
                args, _ = parser.parse_known_args(argv)

                for token in sys.argv[1:]:
                    if token == "--pubchem-enable" or token.startswith("--pubchem-enable="):
                        raise SystemExit("unexpected --pubchem-enable flag present")

                if args.pubchem_enable is not False:
                    raise SystemExit(
                        f"expected pubchem_enable to be False, got {args.pubchem_enable!r}"
                    )

                result_path = os.environ.get("PUBCHEM_TEST_RESULT_PATH")
                if result_path:
                    Path(result_path).write_text(
                        json.dumps(
                            {
                                "pubchem_enable": args.pubchem_enable,
                                "env": os.environ.get("CHEMBL_DA_PUBCHEM_ENABLE"),
                            }
                        ),
                        encoding="utf-8",
                    )

                if args.pubchem_enable:
                    _pc_init_session()

                return 0


            if __name__ == "__main__":
                raise SystemExit(main())
            """
        ),
        encoding="utf-8",
    )

    result_path = tmp_path / "pubchem_result.json"
    monkeypatch.setenv("PUBCHEM_TEST_RESULT_PATH", str(result_path))
    monkeypatch.setenv("CHEMBL_DA_PUBCHEM_ENABLE", "0")

    custom_stages = tuple(
        get_data.Stage(stage.name, str(stub_script) if stage.name == "testitem" else stage.script)
        for stage in get_data.STAGES
    )
    monkeypatch.setattr(get_data, "STAGES", custom_stages)
    monkeypatch.setattr(get_data, "LOGS_DIR", tmp_path / "logs")
    output_dir = (tmp_path / "output").resolve()
    monkeypatch.setattr(get_data, "OUTPUT_DIR", output_dir)
    monkeypatch.setattr(get_data, "CANONICAL_OUTPUT_DIR", output_dir)

    argv = [
        "--skip",
        "document",
        "target",
        "assay",
        "activity",
        "--no-pubchem-enable",
    ]

    exit_code = get_data.main(argv)
    assert exit_code == 0
    assert result_path.exists()

    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["pubchem_enable"] is False
    assert payload["env"] == "0"


_ASSAY_DICTIONARY_PATH = (
    Path(__file__).resolve().parents[1]
    / "resources"
    / "pipeline_inputs"
    / "assay_dictionary.csv"
)


@lru_cache(maxsize=1)
def _load_assay_dictionary() -> pd.DataFrame:
    lookup = pd.read_csv(_ASSAY_DICTIONARY_PATH)
    lookup = lookup.copy()
    lookup["assay_chembl_id"] = lookup["assay_chembl_id"].astype("string").str.strip()
    lookup["assay_strain"] = lookup["assay_strain"].astype("string")
    lookup["assay_group"] = lookup["assay_group"].astype("string")
    lookup["accession"] = lookup["accession"].astype("string")
    lookup["year"] = pd.to_numeric(lookup["year"], errors="coerce").astype("Int64")
    return lookup
