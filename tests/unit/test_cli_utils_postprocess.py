"""Tests for automatic postprocessing hooks in :mod:`library.cli.utils`."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from library.cli.pipeline_definition import PipelineDefinition
from library.cli.utils import _PostprocessHandlers, run_pipeline
from library.postprocess.common import PostprocessingPipelineResult
from library.reporting.run_manifest import PipelineOutputReport


class _LoggerStub:
    def __init__(self) -> None:
        self.messages: list[tuple[str, dict[str, object]]] = []

    def info(self, message: str, **kwargs: object) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def warning(
        self, message: str, **kwargs: object
    ) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def error(
        self, message: str, **kwargs: object
    ) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def exception(
        self, message: str, **kwargs: object
    ) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))


def _writer_stub(
    chunks, destination: Path, col_order, key_cols
):  # pragma: no cover - helper
    collected = [frame.copy() for frame in chunks]
    if collected:
        result = pd.concat(collected, ignore_index=True)
    else:
        result = pd.DataFrame()
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, index=False)
    return destination


def _fetcher_stub():  # pragma: no cover - helper
    yield pd.DataFrame({"id": [1]})


def test_run_pipeline__triggers_postprocessing_on_known_table(monkeypatch, tmp_path):
    captured: dict[str, object] = {}

    class _MetricsStub:
        pipeline_version = "stub-version"

        def summary(self) -> dict[str, int]:
            return {"rows": 0, "columns": 0}

    def _fake_postprocess(table: str, input_path: Path, output_path: Path, config):
        captured["table"] = table
        captured["input"] = Path(input_path)
        captured["output"] = Path(output_path)
        output_path = Path(output_path)
        output_path.write_text("postprocessed", encoding="utf-8")
        report_path = output_path.parent / "activities.postprocess.report.json"
        report_path.write_text("{}", encoding="utf-8")
        return PostprocessingPipelineResult(
            dataframe=pd.DataFrame(),
            metrics=_MetricsStub(),
            output_path=output_path,
            report_path=report_path,
        )

    runtime_stub = SimpleNamespace(
        pipeline_config_factory=lambda **kwargs: SimpleNamespace(**kwargs),
        run_pipeline=_fake_postprocess,
        load_pipeline_config=lambda table, override: SimpleNamespace(
            table=table, override=override
        ),
        get_csv_runtime_config=lambda _cfg: SimpleNamespace(
            separator=",", encoding="utf-8", chunksize=1
        ),
        handlers={
            "activities": _PostprocessHandlers(
                runner=lambda *args, **kwargs: (pd.DataFrame(), _MetricsStub()),
                validator=lambda df: df,
                schema=object(),
            )
        },
    )

    monkeypatch.setattr(
        "library.cli.utils._load_postprocess_runtime",
        lambda: runtime_stub,
    )

    def _fake_finalise(**kwargs):
        captured["finalise_kwargs"] = kwargs
        csv_path = Path(kwargs["csv_path"])
        stats = {
            "rows_total": kwargs.get("rows_total", 0),
            "rows_kept": kwargs.get("rows_kept", 0),
            "rows_dropped": kwargs.get("rows_total", 0) - kwargs.get("rows_kept", 0),
            "output_sha256": "stub",
        }
        extra_stats = kwargs.get("stats_extra") or {}
        stats.update(extra_stats)
        return PipelineOutputReport(
            csv_path=csv_path,
            stats=stats,
            meta_path=csv_path.with_name(csv_path.name + ".meta.yaml"),
            meta_sha256="meta",
        )

    monkeypatch.setattr("library.cli.utils.finalise_csv_output", _fake_finalise)

    definition = PipelineDefinition(
        schema=None,
        schema_name="DummySchema",
        writer=_writer_stub,
        validators=(),
        metadata_hooks=(),
        command="test",
        config_snapshot={},
        inputs={},
        key_columns=(),
        table_quality=lambda path: None,
    )

    output_path = tmp_path / "output.activity_sample.csv"
    failure_path = tmp_path / "failures.csv"

    exit_code = run_pipeline(
        definition=definition,
        fetcher=_fetcher_stub,
        output_path=output_path,
        failure_path=failure_path,
        cfg=None,
        logger=_LoggerStub(),
        postprocess_enabled=True,
    )

    assert exit_code == 0
    assert captured["table"] == "activities"
    assert captured["input"] == output_path
    assert captured["output"].name == "output_postprocessed.activity_sample.csv"
    assert captured["output"].exists()
    finalise_kwargs = captured["finalise_kwargs"]
    assert finalise_kwargs is not None
    extra_metadata = finalise_kwargs.get("extra_metadata") or {}
    assert extra_metadata.get("output_postprocessed") == str(captured["output"])
    assert "postprocess_metrics" in extra_metadata


def test_run_pipeline__skips_postprocessing_when_disabled(monkeypatch, tmp_path):
    def _unexpected_runtime():
        raise AssertionError("postprocess runtime should not load when disabled")

    monkeypatch.setattr(
        "library.cli.utils._load_postprocess_runtime",
        _unexpected_runtime,
    )

    definition = PipelineDefinition(
        schema=None,
        schema_name="DummySchema",
        writer=_writer_stub,
        validators=(),
        metadata_hooks=(),
        command="test",
        config_snapshot={},
        inputs={},
        key_columns=(),
        table_quality=lambda path: None,
    )

    output_path = tmp_path / "output.activity_sample.csv"
    failure_path = tmp_path / "failures.csv"

    captured: dict[str, object] = {}

    def _fake_finalise(**kwargs):
        captured["finalise_kwargs"] = kwargs
        csv_path = Path(kwargs["csv_path"])
        stats = {
            "rows_total": kwargs.get("rows_total", 0),
            "rows_kept": kwargs.get("rows_kept", 0),
            "rows_dropped": kwargs.get("rows_total", 0) - kwargs.get("rows_kept", 0),
            "output_sha256": "stub",
        }
        return PipelineOutputReport(
            csv_path=csv_path,
            stats=stats,
            meta_path=csv_path.with_name(csv_path.name + ".meta.yaml"),
            meta_sha256="meta",
        )

    monkeypatch.setattr("library.cli.utils.finalise_csv_output", _fake_finalise)

    logger = _LoggerStub()
    exit_code = run_pipeline(
        definition=definition,
        fetcher=_fetcher_stub,
        output_path=output_path,
        failure_path=failure_path,
        cfg=None,
        logger=logger,
        postprocess_enabled=False,
    )

    assert exit_code == 0
    assert any(
        message == "[INFO] Postprocessing skipped (flag --postprocess not set)"
        for message, _ in logger.messages
    )
    finalise_kwargs = captured.get("finalise_kwargs") or {}
    extra_metadata = finalise_kwargs.get("extra_metadata") or {}
    assert "output_postprocessed" not in extra_metadata
