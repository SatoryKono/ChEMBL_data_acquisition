"""CLI tests for :mod:`library.utils.cli_tools.pipeline_targets_main`."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import json
import pandas as pd
import pytest
from pandera.errors import SchemaErrors

from library.config import Config
from library.pipelines.target.pipeline import PipelineResult
from library.utils.cli_tools import pipeline_targets_main as cli
from library.schemas import TargetsSchema


class _DummyLogger:
    def __init__(
        self,
        context: dict[str, Any] | None = None,
        storage: list[tuple[str, dict[str, Any]]] | None = None,
    ) -> None:
        self._context: dict[str, Any] = context or {}
        self._records: list[tuple[str, dict[str, Any]]] = (
            storage if storage is not None else []
        )

    def bind(self, **ctx: Any) -> _DummyLogger:  # pragma: no cover - trivial
        merged = {**self._context, **ctx}
        return _DummyLogger(merged, self._records)

    def info(
        self, event: str, *args: Any, **kwargs: Any
    ) -> None:  # pragma: no cover - trivial
        record = {**self._context, **kwargs}
        self._records.append((event, record))

    def error(
        self, event: str, *args: Any, **kwargs: Any
    ) -> None:  # pragma: no cover - trivial
        record = {**self._context, **kwargs}
        self._records.append((event, record))

    @property
    def records(self) -> list[tuple[str, dict[str, Any]]]:  # pragma: no cover - trivial
        return list(self._records)


def test_cli_forwards_batch_size(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **_: Any,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--batch-size",
        "25",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 25
    assert captured["chunks"] == [["CHEMBL1"]]
    assert captured["written_path"] == output_csv
    assert list(captured["written_df"]["target_chembl_id"]) == ["CHEMBL1"]
    assert any(
        event == "pipeline_start" and rec.get("stage") == "pipeline"
        for event, rec in dummy_logger.records
    )
    assert any(
        event == "pipeline_done"
        and rec.get("stage") == "pipeline"
        and rec.get("exit_code") == 0
        for event, rec in dummy_logger.records
    )


def test_cli_limit_restricts_rows(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\nCHEMBL3\n", encoding="utf8"
    )
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **_: Any,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--limit",
        "2",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 100
    assert captured["chunks"] == [["CHEMBL1", "CHEMBL2"]]
    assert captured["written_path"] == output_csv


def test_cli_limit_allows_zero(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text(
        "target_chembl_id\nCHEMBL1\nCHEMBL2\n", encoding="utf8"
    )
    output_csv = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        captured["chunks"] = [list(chunk) for chunk in chunk_iterator()]
        captured["batch_size"] = batch_size
        return PipelineResult(
            chembl=pd.DataFrame({"target_chembl_id": pd.Series(dtype="string")})
        )

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **_: Any,
    ) -> Path:
        captured["written_path"] = Path(path)
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["written_chunks"] = chunks
        captured["written_df"] = (
            pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        )
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),

        "--limit",
        "0",
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    assert captured["batch_size"] == 100
    assert captured["chunks"] == []
    assert captured["written_path"] == output_csv
    assert captured["written_df"].empty


def test_cli_parses_output_arguments(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    captured: dict[str, cli.PipelineConfig] = {}

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(cli, "print_config", lambda cfg: None)

    def fake_run(cfg: Config, options: cli.PipelineConfig) -> int:
        captured["options"] = options
        return 0

    monkeypatch.setattr(cli, "run", fake_run)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
        "--raw-format",
        "parquet",
        "--id-cols",
        "target_chembl_id",
        "pref_name",
        "--no-reindex-raw",
        "--no-normalize-at-export",
    ]
    exit_code = cli.main(args)

    assert exit_code == 0
    assert "options" in captured
    options = captured["options"]
    assert options.raw_format == "parquet"
    assert options.id_cols == ["target_chembl_id", "pref_name"]
    assert options.no_reindex_raw is True
    assert options.normalize_at_export is False



def test_cli_does_not_print_config_when_flag_missing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    output_csv = tmp_path / "out.csv"

    def fake_run_pipeline(
        chunk_iterator: Any,
        cfg: Config,
        *,
        chembl_fetcher: Any,
        batch_size: int | None,
        **_: Any,
    ) -> PipelineResult:
        next(chunk_iterator())
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    def fake_write_csv(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path | str,
        *,
        cfg: Config,
        sep: str | None = None,
        encoding: str | None = None,
        **_: Any,
    ) -> Path:
        if not isinstance(data, pd.DataFrame):
            list(data)
        return Path(path)

    dummy_logger = _DummyLogger()
    monkeypatch.setattr(cli, "configure_logger", lambda cfg: dummy_logger)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a: Config())
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)

    def fail_print_config(_: Config) -> None:
        raise AssertionError("print_config should not be called")

    monkeypatch.setattr(cli, "print_config", fail_print_config)
    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "write_csv", fake_write_csv)

    args = [
        "--input",
        str(input_csv),
        "--output",
        str(output_csv),
    ]
    exit_code = cli.main(args)
    assert exit_code == 0
    captured = capsys.readouterr()
    assert captured.out == ""


def test_cached_chembl_fetch_streams_chunks(monkeypatch: pytest.MonkeyPatch) -> None:
    chunks = iter(
        [
            ["CHEMBL1", "CHEMBL2"],
            ["CHEMBL3"],
        ]
    )

    stream = cli._cached_chembl_fetch(chunks, Config())
    assert hasattr(stream, "__iter__")

    frames = list(stream)
    assert len(frames) == 2
    assert [list(frame["target_chembl_id"]) for frame in frames] == [
        ["CHEMBL1", "CHEMBL2"],
        ["CHEMBL3"],
    ]
    assert all((frame["source"] == "chembl").all() for frame in frames)


# ---------------------------------------------------------------------------
# New helper and integration tests
# ---------------------------------------------------------------------------

def test_raw_dump_keeps_original_columns_and_order(tmp_path: Path, cfg: Config) -> None:
    frames = [
        pd.DataFrame({"pref_name": ["Alpha"], "target_chembl_id": ["CHEMBL1"]}),
        pd.DataFrame({"pref_name": ["Beta"], "target_chembl_id": ["CHEMBL2"]}),
    ]
    raw_path = tmp_path / "raw.csv"

    cli._raw_dump(frames, raw_path, cfg=cfg, sep=",", encoding="utf8")

    written = raw_path.read_text(encoding="utf8").splitlines()
    assert written[0] == "pref_name,target_chembl_id"
    raw_df = pd.read_csv(raw_path)
    assert list(raw_df.columns) == ["pref_name", "target_chembl_id"]
    assert list(raw_df["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]


def test_id_placeholders_replaced_in_id_cols_only() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["-", "CHEMBL2"],
            "pref_name": ["-", "Name"],
            "uniprot_id_primary": ["-", "P12345"],
        }
    )

    cleaned = cli._replace_id_placeholders(frame)

    assert list(cleaned["target_chembl_id"]) == ["", "CHEMBL2"]
    assert list(cleaned["uniprot_id_primary"]) == ["", "P12345"]
    assert list(cleaned["pref_name"]) == ["-", "Name"]


def test_validate_and_write_enforces_TargetsSchema(
    tmp_path: Path, cfg: Config
) -> None:
    valid = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "pref_name": ["Alpha"],
            "uniprot_id_primary": ["-"],
        }
    )
    invalid = pd.DataFrame({"pref_name": ["Alpha"]})

    final_path = tmp_path / "final.csv"
    cli._validate_and_write(valid, final_path, cfg=cfg)

    written = pd.read_csv(final_path, dtype=str, keep_default_na=False)
    assert written.loc[0, "uniprot_id_primary"] == ""
    TargetsSchema.validate(written)

    with pytest.raises(SchemaErrors):
        cli._validate_and_write(invalid, tmp_path / "broken.csv", cfg=cfg)


def test_backward_compatibility_out_alias(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\n", encoding="utf8")
    final_path = tmp_path / "final.csv"

    def fake_run_pipeline(*_: Any, **__: Any) -> PipelineResult:
        return PipelineResult(chembl=pd.DataFrame({"target_chembl_id": ["CHEMBL1"]}))

    recorded: dict[str, Path] = {}

    def fake_validate(
        df: pd.DataFrame | Iterable[pd.DataFrame], path: Path, **_: Any
    ) -> Path:
        recorded["path"] = path
        if not isinstance(df, pd.DataFrame):
            list(df)
        return path

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "_validate_and_write", fake_validate)
    monkeypatch.setattr(cli, "_raw_dump", lambda *a, **k: None)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a, **k: cfg)
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)

    exit_code = cli.main(["--input", str(input_csv), "--out", str(final_path)])

    assert exit_code == 0
    assert recorded["path"] == final_path.resolve()


def test_end_to_end_cli_raw_and_final_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    fixture_path = (
        Path(__file__).resolve().parent / "fixtures" / "chembl_targets_response.json"
    )
    records = json.loads(fixture_path.read_text(encoding="utf8"))
    columns = list(records[0].keys())
    fixture_df = pd.DataFrame(records)[columns]

    def fake_run_pipeline(*_: Any, **__: Any) -> PipelineResult:
        return PipelineResult(chembl=fixture_df.copy())

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli.cli, "apply_config_overrides", lambda *a, **k: cfg)
    monkeypatch.setattr(cli, "ensure_dirs", lambda cfg: None)

    input_csv = tmp_path / "targets.csv"
    input_csv.write_text("target_chembl_id\nCHEMBL1\nCHEMBL2\n", encoding="utf8")
    raw_out = tmp_path / "raw.csv"
    final_out = tmp_path / "final.csv"

    exit_code = cli.main(
        [
            "--input",
            str(input_csv),
            "--raw-out",
            str(raw_out),
            "--final-out",
            str(final_out),
        ]
    )

    assert exit_code == 0

    raw_df = pd.read_csv(raw_out)
    assert list(raw_df.columns) == columns
    assert list(raw_df["target_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]

    final_df = pd.read_csv(final_out, dtype=str, keep_default_na=False)
    TargetsSchema.validate(final_df)
    assert final_df.loc[0, "uniprot_id_primary"] == ""


def test_run_threads_output_flags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    raw_out = tmp_path / "raw.csv"
    final_out = tmp_path / "final.csv"

    base_columns = ["target_chembl_id", "pref_name"]
    metadata_columns = list(cli._PIPELINE_METADATA_COLUMNS)
    payload = {
        "target_chembl_id": ["CHEMBL1"],
        "pref_name": ["Name"],
    }
    for column in metadata_columns:
        payload[column] = [f"meta_{column}"]
    frame = pd.DataFrame(payload)

    def fake_run_pipeline(*_: Any, **__: Any) -> PipelineResult:
        return PipelineResult(chembl=[frame.copy()])

    captured: dict[str, Any] = {"raw_chunks": []}

    class DummyWriter:
        def __init__(
            self,
            path: Path,
            *,
            raw_format: str,
            reindex_columns: bool,
            **kwargs: Any,
        ) -> None:
            captured["writer_args"] = {
                "path": path,
                "raw_format": raw_format,
                "reindex": reindex_columns,
            }
            self.path = path

        def write(self, chunk: pd.DataFrame) -> None:
            captured["raw_chunks"].append(chunk.copy())

        def close(self) -> Path:
            captured["writer_closed"] = True
            return self.path

    def fake_validate(
        data: pd.DataFrame | Iterable[pd.DataFrame],
        path: Path,
        *,
        id_cols: Iterable[str] | None,
        normalize: bool,
        **kwargs: Any,
    ) -> Path:
        captured["validate_path"] = path
        captured["validate_id_cols"] = list(id_cols or [])
        captured["validate_normalize"] = normalize
        chunks: list[pd.DataFrame]
        if isinstance(data, pd.DataFrame):
            chunks = [data.copy()]
        else:
            chunks = [chunk.copy() for chunk in data]
        captured["final_chunks"] = chunks
        return path

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "_RawStreamWriter", DummyWriter)
    monkeypatch.setattr(cli, "_validate_and_write", fake_validate)

    options = cli.PipelineConfig(
        input_csv=tmp_path / "targets.csv",
        final_out=final_out,
        raw_out=raw_out,
        column="target_chembl_id",
        raw_format="csv",
        id_cols=["target_chembl_id", "pref_name"],
        no_reindex_raw=True,
        normalize_at_export=False,
    )

    exit_code = cli.run(cfg, options)
    assert exit_code == 0
    assert captured["writer_args"]["raw_format"] == "csv"
    assert captured["writer_args"]["reindex"] is False
    assert captured["writer_closed"] is True
    assert captured["validate_path"] == final_out
    assert captured["validate_id_cols"] == ["target_chembl_id", "pref_name"]
    assert captured["validate_normalize"] is False
    assert captured["raw_chunks"]
    raw_columns = [list(chunk.columns) for chunk in captured["raw_chunks"]]
    assert raw_columns == [base_columns for _ in raw_columns]
    assert captured["final_chunks"]
    final_columns = [list(chunk.columns) for chunk in captured["final_chunks"]]
    assert final_columns == [base_columns for _ in final_columns]


def test_run_skips_validate_for_parquet_raw_final_match(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    destination = tmp_path / "targets.parquet"
    frame = pd.DataFrame({"target_chembl_id": ["CHEMBL1"]})

    def fake_run_pipeline(*_: Any, **__: Any) -> PipelineResult:
        return PipelineResult(chembl=[frame.copy()])

    class DummyWriter:
        def __init__(self, path: Path, **_: Any) -> None:
            self.path = path
            self.written: list[pd.DataFrame] = []

        def write(self, chunk: pd.DataFrame) -> None:
            self.written.append(chunk.copy())

        def close(self) -> Path:
            return self.path

    monkeypatch.setattr(cli, "run_pipeline", fake_run_pipeline)
    monkeypatch.setattr(cli, "_RawStreamWriter", DummyWriter)

    def fail_validate(*_: Any, **__: Any) -> Path:
        raise AssertionError("_validate_and_write should not be invoked")

    monkeypatch.setattr(cli, "_validate_and_write", fail_validate)

    options = cli.PipelineConfig(
        input_csv=tmp_path / "targets.csv",
        final_out=destination,
        raw_out=destination,
        column="target_chembl_id",
        raw_format="parquet",
        normalize_at_export=False,
    )

    exit_code = cli.run(cfg, options)
    assert exit_code == 0


def test_streaming_pipeline_keeps_chunk_bound(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Config
) -> None:
    total = 1_000
    chunk_size = 25
    input_csv = tmp_path / "targets.csv"
    raw_out = tmp_path / "raw.csv"
    final_out = tmp_path / "final.csv"

    with input_csv.open("w", encoding="utf8") as handle:
        handle.write("target_chembl_id\n")
        for index in range(total):
            handle.write(f"CHEMBL{index}\n")

    recorded_final_sizes: list[int] = []
    recorded_raw_sizes: list[int] = []

    import library.io.writers as io_writers

    original_chunks = io_writers.write_csv_chunks_deterministic

    def tracking_write_chunks(
        chunks: Iterable[pd.DataFrame], path: Path, *args: Any, **kwargs: Any
    ) -> Path:
        def _tracking() -> Iterator[pd.DataFrame]:
            for chunk in chunks:
                recorded_final_sizes.append(len(chunk))
                yield chunk

        return original_chunks(_tracking(), path, *args, **kwargs)

    monkeypatch.setattr(
        io_writers,
        "write_csv_chunks_deterministic",
        tracking_write_chunks,
    )

    original_raw_write = cli._RawStreamWriter.write

    def tracking_raw_write(self: cli._RawStreamWriter, chunk: pd.DataFrame) -> None:
        recorded_raw_sizes.append(len(chunk))
        original_raw_write(self, chunk)

    monkeypatch.setattr(cli._RawStreamWriter, "write", tracking_raw_write)

    options = cli.PipelineConfig(
        input_csv=input_csv,
        final_out=final_out,
        raw_out=raw_out,
        chunk_size=chunk_size,
        batch_size=chunk_size,
        column="target_chembl_id",
        encoding="utf8",
        sep=",",
    )

    exit_code = cli.run(cfg, options)
    assert exit_code == 0

    assert recorded_final_sizes
    assert max(recorded_final_sizes) <= chunk_size
    assert recorded_raw_sizes
    assert max(recorded_raw_sizes) <= chunk_size
