"""Tests for automatic postprocessing hooks in :mod:`library.cli.utils`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from library.cli.pipeline_definition import PipelineDefinition
from library.cli.utils import run_pipeline
from library.postprocess.common import PostprocessingPipelineResult


class _LoggerStub:
    def __init__(self) -> None:
        self.messages: list[tuple[str, dict[str, object]]] = []

    def info(self, message: str, **kwargs: object) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def warning(self, message: str, **kwargs: object) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def error(self, message: str, **kwargs: object) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))

    def exception(self, message: str, **kwargs: object) -> None:  # pragma: no cover - helper
        self.messages.append((message, kwargs))


def _writer_stub(chunks, destination: Path, col_order, key_cols):  # pragma: no cover - helper
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
    captured: dict[str, Path] = {}

    def _fake_postprocess(table: str, input_path: Path, output_path: Path, config):
        captured["table"] = table
        captured["input"] = Path(input_path)
        captured["output"] = Path(output_path)
        output_path = Path(output_path)
        output_path.write_text("postprocessed", encoding="utf-8")
        return PostprocessingPipelineResult(
            dataframe=pd.DataFrame(),
            metrics=None,
            output_path=output_path,
        )

    monkeypatch.setattr(
        "library.cli.utils.run_postprocessing_pipeline",
        _fake_postprocess,
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

    output_path = tmp_path / "output.activities_sample.csv"
    failure_path = tmp_path / "failures.csv"

    exit_code = run_pipeline(
        definition=definition,
        fetcher=_fetcher_stub,
        output_path=output_path,
        failure_path=failure_path,
        cfg=None,
        logger=_LoggerStub(),
    )

    assert exit_code == 0
    assert captured["table"] == "activities"
    assert captured["input"] == output_path
    assert captured["output"].name == "output_postprocessed.activities_sample.csv"
    assert captured["output"].exists()
