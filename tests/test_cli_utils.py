from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd
import pytest

from library.cli_utils import build_parser, run_pipeline
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
        "--config",
        "--sep",
        "--encoding",
        "--col-order",
        "--key-cols",
        "--chunk-size",
        "--merge-chunk-size",
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


def test_run_pipeline_applies_hooks_and_writes(tmp_path: Path) -> None:
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


def test_run_pipeline_writes_failure_cases(tmp_path: Path) -> None:
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
        frames = list(chunks)
        df = (
            pd.concat(frames, ignore_index=True)
            if frames
            else pd.DataFrame(columns=col_order)
        )
        df.to_csv(destination, index=False)
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
    )

    assert exit_code == 1
    assert failure_path.exists()
    failure_csv = failure_path.read_text()
    assert "column" in failure_csv
    assert "value" in failure_csv


def test_run_pipeline_missing_required_columns(tmp_path: Path) -> None:
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
    )

    assert exit_code == 1
    assert not output.exists()
    assert not failure_path.exists()
