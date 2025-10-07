"""Unit tests for the document pipeline programmatic entrypoint."""

from __future__ import annotations

import copy
import sys
from types import ModuleType, SimpleNamespace
from pathlib import Path

import pytest

from library.pipelines.document import DocumentPipelineOptions, run_pipeline


class _DummySection:
    def __init__(self, *, limit: int | None = None, offset: int = 0) -> None:
        self.limit = limit
        self.offset = offset

    def model_copy(self, update: dict[str, object] | None = None) -> "_DummySection":
        update = update or {}
        limit_value = update.get("limit", self.limit)
        if not isinstance(limit_value, int) and limit_value is not None:
            limit_value = self.limit
        offset_value = update.get("offset", self.offset)
        if not isinstance(offset_value, int):
            offset_value = self.offset
        return _DummySection(limit=limit_value, offset=offset_value)


class _DummyConfig:
    def __init__(self) -> None:
        document = SimpleNamespace(
            chembl=_DummySection(),
            pubmed=_DummySection(),
            all=_DummySection(),
        )
        self.sources = SimpleNamespace(
            chembl=SimpleNamespace(pipelines=SimpleNamespace(document=document))
        )

    def model_copy(self, *, deep: bool = False) -> "_DummyConfig":
        # The real config uses pydantic's copy mechanism. A deepcopy is sufficient
        # for these tests and keeps the API compatible.
        return copy.deepcopy(self)


def _install_fake_cli(monkeypatch: pytest.MonkeyPatch, exit_code: int) -> None:
    module = ModuleType("scripts.get_document_data")

    def _runner(*_: object, **__: object) -> int:
        return exit_code

    setattr(module, "run_all", _runner)
    setattr(module, "run_chembl", _runner)
    setattr(module, "run_pubmed", _runner)
    monkeypatch.setitem(sys.modules, "scripts.get_document_data", module)


@pytest.mark.unit
def test_run_pipeline__sets_written_flag_on_success(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _install_fake_cli(monkeypatch, exit_code=0)
    cfg = _DummyConfig()
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id,name,smiles\n", encoding="utf-8")
    options = DocumentPipelineOptions(
        input_csv=input_csv,
        output_csv=tmp_path / "output.csv",
        mode="all",
    )

    result = run_pipeline(cfg, options)

    assert result.exit_code == 0
    assert result.written is True


@pytest.mark.unit
def test_run_pipeline__retains_written_flag_on_skip(tmp_path: Path) -> None:
    cfg = _DummyConfig()
    output_path = tmp_path / "output.csv"
    output_path.write_text("existing", encoding="utf-8")
    options = DocumentPipelineOptions(
        input_csv=tmp_path / "input.csv",
        output_csv=output_path,
        mode="all",
        skip_existing=True,
    )

    result = run_pipeline(cfg, options)

    assert result.executed is False
    assert result.written is False
