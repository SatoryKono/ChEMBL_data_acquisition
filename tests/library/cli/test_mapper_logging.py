import argparse
import importlib
import io
import json
import logging
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

import pandas as pd

from library.cli import LoggerConfig, configure_logger
from library.config import Config
from library.utils.cli_tools import mapper_main


def test_mapper_library_has_no_logging_side_effect() -> None:
    root = logging.getLogger()
    for handler in list(root.handlers):
        root.removeHandler(handler)
    assert not root.handlers
    module_name = "library.integration.mapper_library"
    if module_name in sys.modules:
        del sys.modules[module_name]
    importlib.import_module(module_name)
    assert not root.handlers


def test_mapper_main_logs_mapping(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    """Mapper CLI emits structured mapping log records."""

    def fake_map(
        ids: list[str],
        cfg: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        assert ids == ["CHEMBL1"]
        assert batch_size == 1
        assert rps == 1.0
        assert max_workers == 1
        return {"CHEMBL1": "P12345"}

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    df = pd.DataFrame({"chembl_id": ["CHEMBL1"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=1,
        rps=1.0,
        workers=1,
    )
    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )
    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0
    records = [json.loads(line) for line in buffer.getvalue().splitlines()]
    summary = [r for r in records if r["event"] == "mapper_summary"]
    assert summary
    payload = summary[0]
    assert payload["total"] == 1
    assert payload["mapped"] == 1
    assert payload["missing"] == 0
    assert "sample_missing" not in payload
    assert not [r for r in records if r["event"] in {"mapped", "uniprot_id_missing"}]


def test_mapper_main_concurrent_preserves_order(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    """Concurrent path retains mapping order and logging semantics."""

    df = pd.DataFrame({"chembl_id": ["CHEMBL1", None, "CHEMBL2", " ", "CHEMBL3"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        captured["ids"] = list(ids)
        captured["batch_size"] = batch_size
        captured["rps"] = rps
        captured["max_workers"] = max_workers
        return {
            "CHEMBL2": "P222",
            "CHEMBL1": "P111",
            "CHEMBL3": None,
        }

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    written: dict[str, Any] = {}

    def fake_write_csv(
        df_out: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str,
        encoding: str,
        key_cols: Any,
    ) -> Path:
        written["df"] = df_out.copy()
        written["path"] = path
        written["key_cols"] = key_cols
        written["chembl_ids"] = df_out["chembl_id"].tolist()
        return path

    monkeypatch.setattr(mapper_main.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=output_path,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=2,
        rps=5.0,
        workers=3,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="DEBUG"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0

    assert captured["ids"] == ["CHEMBL1", "CHEMBL2", "CHEMBL3"]
    assert captured["batch_size"] == 2
    assert captured["rps"] == 5.0
    assert captured["max_workers"] == 3

    output_df = written["df"]
    chembl_ids = [
        value if isinstance(value, str) and value.strip() else None
        for value in written["chembl_ids"]
    ]
    assert chembl_ids == ["CHEMBL1", None, "CHEMBL2", "CHEMBL3"]
    mapped_values = [
        value if isinstance(value, str) and value else None
        for value in output_df["mapping_uniprot_id"].tolist()
    ]
    assert mapped_values == ["P111", None, "P222", None]
    assert written["path"] == output_path
    assert written["key_cols"] is None

    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    status_events = [
        (record["event"], record.get("chembl_id"))
        for record in records
        if record["event"] in {"mapped", "uniprot_id_missing"}
    ]
    assert status_events == [
        ("mapped", "CHEMBL1"),
        ("mapped", "CHEMBL2"),
        ("uniprot_id_missing", "CHEMBL3"),
    ]
    summary = [r for r in records if r["event"] == "mapper_summary"]
    assert summary
    payload = summary[0]
    assert payload["total"] == 3
    assert payload["mapped"] == 2
    assert payload["missing"] == 1
    assert payload["sample_missing"] == ["CHEMBL3"]


def test_mapper_main_deduplicates_ids_before_mapping(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    df = pd.DataFrame({"chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    calls: dict[str, Any] = {"count": 0, "ids": None}

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        calls["count"] += 1
        calls["ids"] = list(ids)
        return {"CHEMBL1": "P111", "CHEMBL2": "P222"}

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    written: dict[str, Any] = {}

    def fake_write_csv(
        df_out: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str,
        encoding: str,
        key_cols: Any,
    ) -> Path:
        written["df"] = df_out.copy()
        written["path"] = path
        written["key_cols"] = key_cols
        return path

    monkeypatch.setattr(mapper_main.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=output_path,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=2,
        rps=5.0,
        workers=3,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="DEBUG"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0

    assert calls["count"] == 1
    assert calls["ids"] == ["CHEMBL1", "CHEMBL2"]

    output_df = written["df"]
    assert output_df["chembl_id"].tolist() == ["CHEMBL1", "CHEMBL1", "CHEMBL2"]
    assert output_df["mapping_uniprot_id"].tolist() == ["P111", "P111", "P222"]
    assert written["path"] == output_path
    assert written["key_cols"] is None

    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    summary = [r for r in records if r["event"] == "mapper_summary"]
    assert summary
    payload = summary[0]
    assert payload["total"] == 3
    assert payload["mapped"] == 3
    assert payload["missing"] == 0
    assert "sample_missing" not in payload


def test_mapper_main_log_each_emits_per_id_logs(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    df = pd.DataFrame({"chembl_id": ["CHEMBL1", "CHEMBL2"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        return {"CHEMBL1": "P111", "CHEMBL2": None}

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,

        chunk_size=1,
        rps=1.0,
        workers=1,
        log_each=True,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0

    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    mapped = [r for r in records if r["event"] == "mapped"]
    assert mapped and mapped[0]["level"] == "INFO"
    missing = [r for r in records if r["event"] == "uniprot_id_missing"]
    assert missing and missing[0]["level"] in {"WARN", "WARNING"}
    summary = [r for r in records if r["event"] == "mapper_summary"]
    assert summary
    payload = summary[0]
    assert payload["total"] == 2
    assert payload["mapped"] == 1
    assert payload["missing"] == 1
    assert payload["sample_missing"] == ["CHEMBL2"]


def test_mapper_main_timeout_logs_missing_summary(
    tmp_path: Path, monkeypatch: Any, cfg: Config
) -> None:
    df = pd.DataFrame({"chembl_id": ["CHEMBL1", "CHEMBL2"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)

    def raise_timeout(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        raise TimeoutError("network timeout")

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", raise_timeout)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=tmp_path / "out.csv",
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=2,
        rps=5.0,
        workers=3,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 1

    records = [json.loads(line) for line in buffer.getvalue().splitlines() if line]
    summary = [record for record in records if record["event"] == "mapper_summary"]
    assert summary
    payload = summary[0]
    assert payload["total"] == 2
    assert payload["mapped"] == 0
    assert payload["missing"] == 2
    assert payload["sample_missing"] == ["CHEMBL1", "CHEMBL2"]
    assert not [
        record
        for record in records
        if record["event"] in {"mapped", "uniprot_id_missing"}
    ]


@pytest.mark.parametrize(
    "keep_markers, expected",
    [
        (False, ["CHEMBL1", "CHEMBL3"]),
        (True, ["CHEMBL1", "CUSTOM_NA", "CHEMBL3"]),
    ],
)
def test_mapper_main_respects_na_marker_configuration(
    tmp_path: Path,
    monkeypatch: Any,
    cfg: Config,
    keep_markers: bool,
    expected: list[str],
) -> None:
    cfg.io.na_markers = ("CUSTOM_NA",)
    cfg.io.keep_na_markers = keep_markers

    df = pd.DataFrame({"chembl_id": ["CHEMBL1", "CUSTOM_NA", "CHEMBL3"]})
    input_path = tmp_path / "in.csv"
    df.to_csv(input_path, index=False)
    output_path = tmp_path / "out.csv"

    captured: dict[str, Any] = {}

    def fake_map(
        ids: list[str],
        cfg_mapping: object,
        *,
        batch_size: int,
        rps: float,
        max_workers: int | None,
    ) -> dict[str, str | None]:
        captured["ids"] = list(ids)
        return {chembl_id: None for chembl_id in ids}

    monkeypatch.setattr(mapper_main, "map_chembl_ids_to_uniprot", fake_map)

    def fake_write_csv(
        df_out: pd.DataFrame,
        path: Path,
        *,
        cfg: Config,
        sep: str,
        encoding: str,
        key_cols: Any,
    ) -> Path:
        return path

    monkeypatch.setattr(mapper_main.io, "write_csv", fake_write_csv)

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=output_path,
        column="chembl_id",
        sep=",",
        encoding="utf8",
        key_cols=None,
        chunk_size=1,
        rps=1.0,
        workers=1,
        log_each=False,
    )

    buffer = io.StringIO()
    configure_logger(LoggerConfig(stream=buffer, level="INFO"))

    cfg_ns = SimpleNamespace(
        io=cfg.io,
        uniprot_mapping=cfg.uniprot_mapping,
        to_dict=lambda: {},
    )

    exit_code = mapper_main.run(cfg_ns, args)
    assert exit_code == 0
    assert captured["ids"] == expected
