"""Unit tests for the document pipeline programmatic entrypoint."""

from __future__ import annotations

import copy
import argparse
from pathlib import Path
from types import SimpleNamespace

import pytest

from library.pipelines.document import DocumentPipelineOptions, run_pipeline
from library.pipelines.document import service as document_service
from library.config.models import ApiCfg, RetryCfg, ChemblCacheCfg
from library.cli.commands import get_document_data as document_cli


class _DummySection:
    def __init__(
        self,
        *,
        limit: int | None = None,
        offset: int = 0,
        **extra: object,
    ) -> None:
        self.limit = limit
        self.offset = offset
        for key, value in extra.items():
            setattr(self, key, value)

    def model_copy(self, update: dict[str, object] | None = None) -> _DummySection:
        update = update or {}
        data = dict(self.__dict__)

        limit_value = update.get("limit", data.get("limit"))
        if not isinstance(limit_value, int) and limit_value is not None:
            limit_value = data.get("limit")
        data["limit"] = limit_value

        offset_value = update.get("offset", data.get("offset", 0))
        if not isinstance(offset_value, int):
            offset_value = data.get("offset", 0)
        data["offset"] = offset_value

        for key, value in update.items():
            if key in {"limit", "offset"}:
                continue
            data[key] = value

        return _DummySection(**data)


class _DummyConfig:
    def __init__(self) -> None:
        document_all = _DummySection(
            column="molecule_chembl_id",
            chunk_size=20,
            sleep=5.0,
            workers=1,
            batch_size=50,
            timeout=90.0,
        )
        document_chembl = _DummySection(
            column="molecule_chembl_id",
            chunk_size=20,
            timeout=90.0,
        )
        document_pubmed = _DummySection(
            column="PMID",
            sleep=5.0,
            workers=1,
            batch_size=100,
        )
        document = SimpleNamespace(
            chembl=document_chembl,
            pubmed=document_pubmed,
            all=document_all,
        )
        self.document = document
        self.sources = SimpleNamespace(
            chembl=SimpleNamespace(pipelines=SimpleNamespace(document=document))
        )
        self.io = SimpleNamespace(
            csv_sep=",",
            csv_encoding="utf-8-sig",
            output_dir=Path.cwd(),
            na_markers=("#N/A",),
            keep_na_markers=False,
        )
        self.local = SimpleNamespace(io=self.io)
        self.api = ApiCfg()
        self.retry = RetryCfg()
        self.chembl = ChemblCacheCfg()

    def model_copy(self, *, deep: bool = False) -> _DummyConfig:
        # The real config uses pydantic's copy mechanism. A deepcopy is sufficient
        # for these tests and keeps the API compatible.
        return copy.deepcopy(self)


def _install_fake_cli(monkeypatch: pytest.MonkeyPatch, exit_code: int) -> None:
    def _runner(cfg: object, args: SimpleNamespace, **_: object) -> int:
        final_out = Path(args.final_out)
        final_out.parent.mkdir(parents=True, exist_ok=True)
        final_out.write_text("identifier\n", encoding="utf-8")

        date_tag = getattr(args, "_standard_date_tag", "") or "19700101"
        table_name = Path(args.final_out).stem or "documents"
        output_dir = Path(getattr(getattr(cfg, "io", SimpleNamespace()), "output_dir", final_out.parent))
        canonical = output_dir / f"output.{table_name}_{date_tag}.csv"
        canonical.parent.mkdir(parents=True, exist_ok=True)
        canonical.write_text("identifier\n", encoding="utf-8")
        return exit_code

    monkeypatch.setattr(
        document_cli,
        "MODE_HANDLERS",
        {"all": _runner, "chembl": _runner, "pubmed": _runner},
        raising=False,
    )
    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None, raising=False)


@pytest.mark.unit
def test_run_pipeline__sets_written_flag_on_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _install_fake_cli(monkeypatch, exit_code=0)
    cfg = _DummyConfig()
    cfg.io.output_dir = tmp_path
    cfg.local.io.output_dir = tmp_path
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


@pytest.mark.unit
def test_run_pipeline__passes_mode_to_cli(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}

    def _runner(cfg: object, args: SimpleNamespace, **_: object) -> int:
        captured["cfg"] = cfg
        captured["args"] = args
        return 0
    monkeypatch.setattr(
        document_cli,
        "MODE_HANDLERS",
        {"all": _runner, "chembl": _runner, "pubmed": _runner},
        raising=False,
    )
    monkeypatch.setattr(document_service, "_MODE_HANDLERS_CACHE", None, raising=False)

    cfg = _DummyConfig()
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("molecule_chembl_id,name,smiles\n", encoding="utf-8")
    options = DocumentPipelineOptions(
        input_csv=input_csv,
        output_csv=tmp_path / "output.csv",
        mode="chembl",
    )

    result = run_pipeline(cfg, options)

    assert result.exit_code == 0
    assert captured
    args = captured["args"]
    assert isinstance(args, argparse.Namespace)
    assert args.mode == "chembl"
    assert args.command == "chembl"
