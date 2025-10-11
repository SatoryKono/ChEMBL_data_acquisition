from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from library.cli.pipeline_definition import PipelineDefinition
from library.cli_utils import RunPipelineResult, run_pipeline
from library.config import IoCfg


def _make_cfg(tmp_path: Path) -> SimpleNamespace:
    io_cfg = IoCfg(
        output_dir=tmp_path,
        exist_ok=True,
        csv_sep=",",
        csv_encoding="utf-8",
        default_date_prefix="20240101",
    )
    system_cfg = SimpleNamespace(doc_quality=SimpleNamespace(fatal_on_error=False))
    return SimpleNamespace(io=io_cfg, system=system_cfg)


@pytest.mark.unit
def test_run_pipeline_result__exposes_attributes() -> None:
    dataset = Path("/tmp/dataset.csv")
    result = RunPipelineResult(0, dataset, None)

    assert int(result) == 0
    assert result == 0
    assert 0 == result
    assert result != 1
    assert not result
    assert hash(result) == hash(int(result))
    assert {result} == {0}
    assert result in {0}
    assert 0 in {result}
    assert result.exit_code == 0
    assert result.dataset_path == dataset
    assert result.artifacts is None


@pytest.mark.unit
def test_run_pipeline__persists_standard_outputs(tmp_path: Path) -> None:
    frame = pd.DataFrame({"identifier": ["row-1", "row-2"], "value": [1, 2]})

    def fetcher() -> list[pd.DataFrame]:
        return [frame]

    writer_called = []

    def writer(
        chunks: list[pd.DataFrame],
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        writer_called.append(True)
        combined = pd.concat(chunks, ignore_index=True)
        combined.to_csv(destination, index=False)
        return destination

    table_quality_called = []

    def table_quality(_: Path) -> None:
        table_quality_called.append(True)

    captured_stats: list[dict[str, object]] = []

    def stats_callback(payload: dict[str, object]) -> None:
        captured_stats.append(dict(payload))

    definition = PipelineDefinition(
        schema=None,
        schema_name="TestSchema",
        writer=writer,
        validators=(),
        metadata_hooks=(),
        command="test",  # command is required by run_pipeline
        config_snapshot={},
        inputs={},
        key_columns=("identifier",),
        table_quality=table_quality,
        stats_callback=stats_callback,
    )

    output_path = tmp_path / "output.documents_20240101.csv"
    failure_path = tmp_path / "failures.csv"
    cfg = _make_cfg(tmp_path)

    result = run_pipeline(
        definition=definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        emit_standard_outputs=True,
        emit_legacy_artifacts=False,
    )

    assert isinstance(result, RunPipelineResult)
    assert int(result) == 0
    assert not writer_called, "legacy writer must not be invoked"
    assert not table_quality_called, "quality hook should be skipped without legacy artefacts"

    artifacts = result.artifacts
    assert artifacts is not None
    assert artifacts.dataset.exists()
    assert artifacts.quality_report.exists()
    assert artifacts.correlation_report.exists()

    dataset = pd.read_csv(artifacts.dataset)
    pd.testing.assert_frame_equal(dataset, frame)

    assert captured_stats, "stats callback must receive payload"
    assert captured_stats[-1]["dataset_path"] == str(artifacts.dataset)


@pytest.mark.unit
def test_run_pipeline__legacy_mode_streams_and_writes_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    frames = [
        pd.DataFrame({"identifier": ["row-1"], "value": [1]}),
        pd.DataFrame({"identifier": ["row-2"], "value": [2]}),
    ]

    def fetcher() -> list[pd.DataFrame]:
        return frames

    writer_calls: list[pd.DataFrame] = []

    def writer(
        chunks: object,
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        assert not isinstance(chunks, list), "chunks must be streamed for legacy mode"
        streamed: list[pd.DataFrame] = []
        for chunk in chunks:  # type: ignore[assignment]
            streamed.append(chunk.copy())
        writer_calls.extend(streamed)
        combined = pd.concat(streamed, ignore_index=True)
        combined.to_csv(destination, index=False)
        return destination

    table_quality_paths: list[Path] = []

    def table_quality(path: Path) -> None:
        table_quality_paths.append(path)

    def forbid_standard_outputs(*_: object, **__: object) -> None:
        raise AssertionError("save_standard_outputs must not be called in legacy-only mode")

    monkeypatch.setattr("library.cli_utils.save_standard_outputs", forbid_standard_outputs)

    definition = PipelineDefinition(
        schema=None,
        schema_name="TestSchema",
        writer=writer,
        validators=(),
        metadata_hooks=(),
        command="test",
        config_snapshot={},
        inputs={},
        key_columns=("identifier",),
        table_quality=table_quality,
    )

    output_path = tmp_path / "output.documents_20240101.csv"
    failure_path = tmp_path / "failures.csv"
    cfg = _make_cfg(tmp_path)

    result = run_pipeline(
        definition=definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        emit_standard_outputs=False,
        emit_legacy_artifacts=True,
    )

    assert isinstance(result, RunPipelineResult)
    assert int(result) == 0
    assert writer_calls, "legacy writer must process streamed chunks"
    assert result.artifacts is None
    assert result.dataset_path == output_path
    assert table_quality_paths == [output_path]

    meta_path = output_path.with_name(output_path.name + ".meta.yaml")
    assert meta_path.exists()

    dataset = pd.read_csv(output_path)
    expected = pd.concat(frames, ignore_index=True)
    pd.testing.assert_frame_equal(dataset, expected)

    assert not failure_path.exists(), "failure artefacts should be cleaned on success"


@pytest.mark.unit
def test_run_pipeline__passes_config_snapshot_to_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    frame = pd.DataFrame({"identifier": ["row-1"], "value": [7]})

    def fetcher() -> list[pd.DataFrame]:
        return [frame]

    def writer(
        chunks: object,
        destination: Path,
        col_order: list[str] | None,
        key_cols: list[str],
    ) -> Path:
        pd.concat(list(chunks), ignore_index=True).to_csv(destination, index=False)
        return destination

    captured_kwargs: list[dict[str, object]] = []

    def fake_write_meta_yaml(csv_path: Path, **kwargs: object) -> Path:
        captured_kwargs.append(dict(kwargs))
        meta_path = Path(csv_path).with_name(Path(csv_path).name + ".meta.yaml")
        meta_path.write_text("{}", encoding="utf-8")
        return meta_path

    monkeypatch.setattr("library.cli_utils.write_meta_yaml", fake_write_meta_yaml)

    config_snapshot = {"io": {"mode": "test"}}

    definition = PipelineDefinition(
        schema=None,
        schema_name="TestSchema",
        writer=writer,
        validators=(),
        metadata_hooks=(),
        command="test",
        config_snapshot=config_snapshot,
        inputs={},
        key_columns=("identifier",),
    )

    output_path = tmp_path / "output.documents_20240101.csv"
    failure_path = tmp_path / "failures.csv"
    cfg = _make_cfg(tmp_path)

    result = run_pipeline(
        definition=definition,
        fetcher=fetcher,
        output_path=output_path,
        failure_path=failure_path,
        cfg=cfg,
        emit_standard_outputs=False,
        emit_legacy_artifacts=True,
    )

    assert int(result) == 0
    assert captured_kwargs, "metadata writer must capture invocation"
    assert captured_kwargs[-1]["config"] == config_snapshot
