from __future__ import annotations

import io
import json
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest
import yaml

from library.config import Config, _mask_secrets, _serialize_paths
from library.cli_utils import build_parser, run_pipeline
from library.logging_setup import LoggerConfig, configure_logger as setup_logger
from schemas import AssaysSchema


def test_cli_utils_flags_and_help() -> None:
    parser = build_parser()
    actions = parser._option_string_actions
    expected = {
        "-h",
        "--help",
        "--log-level",
        "--input",
        "--output",
        "--out",
        "--config",
        "--sep",
        "--encoding",
        "--col-order",
        "--key-cols",
        "--chunk-size",
        "--merge-chunk-size",
        "--base-path",
        "--date",
        "--force",
        "--input-dir",
        "--output-dir",
        "--print-config",
        "--skip-existing",
    }
    assert set(actions) == expected
    assert actions["--log-level"].help == "Logging level"
    assert actions["--input"].help == "Input CSV file"
    assert (
        actions["--output"].help
        == "Destination CSV file (default: output.<stem>_<YYYYMMDD>.csv)"
    )
    assert actions["--sep"].help == "CSV delimiter"
    assert actions["--encoding"].help == "File encoding"
    assert actions["--col-order"].help == "Preferred column order"
    assert actions["--key-cols"].help == "Columns used for sorting"
    assert actions["--chunk-size"].help == "Number of rows read per chunk"
    assert (
        actions["--merge-chunk-size"].help
        == "Rows loaded per temporary file during merge"
    )
    assert parser.description is not None
    assert parser.description.startswith(
        "CLI wrapper for :func:`write_csv_deterministic`"
    )


class _ValidationResult:
    def __init__(self, data: pd.DataFrame, failure_cases: pd.DataFrame) -> None:
        self.data = data
        self.failure_cases = failure_cases


def test_run_pipeline_applies_hooks_and_writes(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"
    frames = [
        pd.DataFrame(
            {
                "assay_chembl_id": ["A1", "A2"],
                "document_chembl_id": ["D1", "D2"],
                "extra": ["x", "y"],
            }
        )
    ]

    def fetcher() -> list[pd.DataFrame]:
        return frames

    hooks = [
        lambda df: df.assign(extra=df["extra"].str.upper()),
        lambda df: df.assign(new_col="value"),
    ]

    def validator(df: pd.DataFrame) -> _ValidationResult:
        return _ValidationResult(df, pd.DataFrame())

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str],
        key_cols: list[str],
    ) -> Path:
        frames = list(chunks)
        if frames:
            df = pd.concat(frames, ignore_index=True)
        else:
            df = pd.DataFrame(columns=col_order)
        df.to_csv(destination, index=False)
        return destination

    quality_called: list[Path] = []

    def quality(path: Path) -> None:
        quality_called.append(path)

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[validator],
        metadata_hooks=hooks,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=quality,
        cfg=cfg,
    )

    assert exit_code == 0
    assert output.exists()
    assert not failure_path.exists()
    assert quality_called == [output]

    written = pd.read_csv(output)
    schema_columns = list(AssaysSchema.columns)
    head = [c for c in schema_columns if c in written.columns]
    tail = sorted(c for c in written.columns if c not in schema_columns)
    assert list(written.columns) == head + tail

    meta_path = output.with_name(output.name + ".meta.yaml")
    assert meta_path.exists()


def test_run_pipeline_removes_outputs_on_table_quality_failure(
    tmp_path: Path, cfg: Config
) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"

    frames = [
        pd.DataFrame(
            {
                "assay_chembl_id": ["A1"],
                "document_chembl_id": ["D1"],
            }
        )
    ]

    def fetcher() -> list[pd.DataFrame]:
        return frames

    def validator(df: pd.DataFrame) -> _ValidationResult:
        return _ValidationResult(df, pd.DataFrame())

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str],
        key_cols: list[str],
    ) -> Path:
        frames = list(chunks)
        pd.concat(frames, ignore_index=True).to_csv(destination, index=False)
        return destination

    quality_calls: list[Path] = []

    def table_quality(path: Path) -> None:
        quality_calls.append(path)
        raise RuntimeError("quality failed")

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[validator],
        metadata_hooks=[lambda df: df],
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=table_quality,
        cfg=cfg,
    )

    meta_path = output.with_name(output.name + ".meta.yaml")

    assert exit_code == 1
    assert quality_calls == [output]
    assert not output.exists()
    assert not meta_path.exists()


def test_run_pipeline_removes_stale_failure_outputs(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"
    failure_meta_path = Path(f"{failure_path}.meta.yaml")

    failure_path.write_text("old", encoding="utf8")
    failure_meta_path.write_text("old", encoding="utf8")

    def fetcher() -> list[pd.DataFrame]:
        return [
            pd.DataFrame(
                {
                    "assay_chembl_id": ["A1"],
                    "document_chembl_id": ["D1"],
                }
            )
        ]

    def validator(df: pd.DataFrame) -> _ValidationResult:
        return _ValidationResult(df, pd.DataFrame())

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str],
        key_cols: list[str],
    ) -> Path:
        frames = list(chunks)
        pd.concat(frames, ignore_index=True).to_csv(destination, index=False)
        return destination

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[validator],
        metadata_hooks=[lambda df: df],
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=lambda path: None,
        cfg=cfg,
    )

    assert exit_code == 0
    assert output.exists()
    assert not failure_path.exists()
    assert not failure_meta_path.exists()


def test_run_pipeline_writes_failure_cases(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"

    def fetcher() -> list[pd.DataFrame]:
        return [
            pd.DataFrame(
                {
                    "assay_chembl_id": ["A1"],
                    "document_chembl_id": ["D1"],
                }
            )
        ]

    def validator(df: pd.DataFrame) -> _ValidationResult:
        failures = pd.DataFrame([{"column": "extra", "value": "bad"}])
        return _ValidationResult(df, failures)

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str],
        key_cols: list[str],
    ) -> Path:
        pytest.fail("writer should not run when validation fails")

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[validator],
        metadata_hooks=[lambda df: df],
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=lambda path: None,
        cfg=cfg,
    )

    assert exit_code == 1
    assert failure_path.exists()
    failure_csv = failure_path.read_text()
    assert "column" in failure_csv
    assert "value" in failure_csv
    assert not output.exists()
    output_meta = output.with_name(output.name + ".meta.yaml")
    assert not output_meta.exists()
    meta_path = Path(str(failure_path) + ".meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf8"))
    expected_config = _mask_secrets(_serialize_paths(cfg.to_dict()))
    normalized_expected = yaml.safe_load(yaml.safe_dump(expected_config))
    assert meta["config"] == normalized_expected


def test_run_pipeline_missing_required_columns(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"

    def fetcher() -> list[pd.DataFrame]:
        return [pd.DataFrame({"document_chembl_id": ["D1"]})]

    hooks = [lambda df: df]

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str],
        key_cols: list[str],
    ) -> Path:
        pytest.fail("writer should not be called when required columns are missing")

    def quality(_: Path) -> None:
        pytest.fail("table quality should not run when required columns are missing")

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[],
        metadata_hooks=hooks,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=quality,
        cfg=cfg,
    )

    assert exit_code == 1
    assert not output.exists()
    assert not failure_path.exists()


def test_run_pipeline_metadata_hook_failure(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"
    frames = [
        pd.DataFrame(
            {
                "assay_chembl_id": ["A1"],
                "document_chembl_id": ["D1"],
            }
        )
    ]

    def fetcher() -> list[pd.DataFrame]:
        return frames

    def hook(_: pd.DataFrame) -> pd.DataFrame:
        raise RuntimeError("hook boom")

    writer_called: list[Path] = []

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        writer_called.append(destination)
        return destination

    buf = io.StringIO()
    logger = setup_logger(LoggerConfig(stream=buf), replace_root=False)

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[],
        metadata_hooks=[hook],
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=lambda path: None,
        cfg=cfg,
        logger=logger,
    )

    assert exit_code == 1
    assert not writer_called
    assert not output.exists()
    assert not failure_path.exists()

    lines = [line for line in buf.getvalue().splitlines() if line.strip()]
    assert lines
    assert all("Traceback" not in line for line in lines)
    records = [json.loads(line) for line in lines]
    hook_errors = [record for record in records if record.get("event") == "metadata_hook_failed"]
    assert hook_errors
    assert hook_errors[-1]["error"] == "hook boom"


def test_run_pipeline_writer_failure(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"

    def fetcher() -> list[pd.DataFrame]:
        return [
            pd.DataFrame(
                {
                    "assay_chembl_id": ["A1"],
                    "document_chembl_id": ["D1"],
                }
            )
        ]

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        raise RuntimeError("writer boom")

    buf = io.StringIO()
    logger = setup_logger(LoggerConfig(stream=buf), replace_root=False)

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=AssaysSchema,
        schema_name="AssaysSchema",
        validators=[],
        metadata_hooks=[lambda df: df],
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=["assay_chembl_id"],
        table_quality=lambda path: None,
        cfg=cfg,
        logger=logger,
    )

    assert exit_code == 1
    assert not output.exists()
    assert not failure_path.exists()

    lines = [line for line in buf.getvalue().splitlines() if line.strip()]
    assert lines
    assert all("Traceback" not in line for line in lines)
    records = [json.loads(line) for line in lines]
    write_errors = [record for record in records if record.get("event") == "write_fail"]
    assert write_errors
    assert write_errors[-1]["error"] == "writer boom"


def test_run_pipeline_table_quality_failure(tmp_path: Path, cfg: Config) -> None:
    output = tmp_path / "assays.csv"
    failure_path = tmp_path / "assays_failure_cases.csv"

    frames = [
        pd.DataFrame(
            {
                "assay_chembl_id": ["A1"],
                "document_chembl_id": ["D1"],
            }
        )
    ]

    def fetcher() -> list[pd.DataFrame]:
        return frames

    def writer(
        chunks: Iterable[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        collected = list(chunks)
        if collected:
            df = pd.concat(collected, ignore_index=True)
        else:
            df = pd.DataFrame(columns=col_order or [])
        df.to_csv(destination, index=False)
        return destination

    buf = io.StringIO()
    logger = setup_logger(LoggerConfig(stream=buf), replace_root=False)

    def failing_quality(_: Path) -> None:
        raise RuntimeError("quality boom")

    exit_code = run_pipeline(
        fetcher=fetcher,
        schema=None,
        schema_name="",
        validators=[],
        metadata_hooks=None,
        writer=writer,
        output_path=output,
        failure_path=failure_path,
        command="pytest",
        config_snapshot={},
        inputs={},
        key_columns=[],
        table_quality=failing_quality,
        cfg=cfg,
        logger=logger,
    )

    assert exit_code == 1
    assert output.exists()
    assert not failure_path.exists()

    lines = [line for line in buf.getvalue().splitlines() if line.strip()]
    assert lines
    records = [json.loads(line) for line in lines]
    quality_errors = [
        record for record in records if record.get("event") == "quality_report_failed"
    ]
    assert quality_errors
    last_record = quality_errors[-1]
    assert last_record["error"] == "quality boom"
    assert last_record.get("traceback")
