from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from pytest import MonkeyPatch

from library.integration import chembl_library as cl
from library import io
from library.config import Config
from scripts import get_activity_data as gad


def _create_config(tmp_path: Path) -> Path:
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "sources:\n"
        "  chembl:\n"
        "    pipelines:\n"
        "      activity:\n"
        "        batch_size: 10\n"
        "    api:\n"
        "      timeout_read: 30\n"
        "local:\n"
        "  io:\n"
        "    csv_sep: '|'\n"
        "    csv_encoding: iso-8859-1\n"
        "  resources:\n"
        "    dictionary_dir: dictionary\n"
        "    iuphar_target_csv: dictionary/_IUPHAR/_IUPHAR_target.csv\n"
        "    iuphar_family_csv: dictionary/_IUPHAR/_IUPHAR_family.csv\n"
        "    uniprot_data_dir: uniprot\n"
        "    targets_type_csv: dictionary/targets_type.csv\n"
        "system:\n"
        "  log:\n"
        "    level: INFO\n"
    )
    return cfg


def _run(
    tmp_path: Path, monkeypatch: MonkeyPatch, extra: list[str]
) -> tuple[int, dict[str, object]]:
    input_csv = tmp_path / "activity.csv"
    input_csv.write_text("activity_id\n1\n")
    config_path = _create_config(tmp_path)
    called: dict[str, object] = {}
    monkeypatch.setattr(io, "read_ids", lambda *a, **k: iter(["1"]))

    def fake_get(
        ids: Sequence[str],
        cfg: Config,
        client: Any,
        chunk_size: int,
        timeout: float,
        **kwargs: object,
    ) -> pd.DataFrame:
        data = list(ids)
        called["batch_size"] = chunk_size
        called["extra_columns"] = kwargs.get("extra_columns")
        return pd.DataFrame({"activity_id": data})

    def fake_write(
        df: pd.DataFrame,
        output: Path,
        *,
        cfg: Config | None = None,
        sep: str | None = None,
        encoding: str | None = None,
        **_: object,
    ) -> Path:
        if sep is None and cfg is not None:
            sep = cfg.io.csv_sep
        if encoding is None and cfg is not None:
            encoding = cfg.io.csv_encoding
        called["sep"] = sep
        called["encoding"] = encoding
        return output

    monkeypatch.setattr(cl, "get_activities", fake_get)
    monkeypatch.setattr(io, "write_csv", fake_write)
    monkeypatch.setattr(
        gad, "analyze_table_quality", lambda df, table_name, **_: None
    )
    monkeypatch.setattr(gad, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gad, "write_meta_yaml", lambda **__: None)
    rc = gad.main(["--config", str(config_path), "--input", str(input_csv), *extra])
    return rc, called


def test_default_config_used(tmp_path: Path, monkeypatch: MonkeyPatch) -> None:
    rc, called = _run(tmp_path, monkeypatch, [])
    assert rc == 0
    assert called["batch_size"] == 10
    assert called["sep"] == "|"
    assert called["encoding"] == "iso-8859-1"


def test_cli_overrides(tmp_path: Path, monkeypatch: MonkeyPatch) -> None:
    rc, called = _run(
        tmp_path,
        monkeypatch,
        ["--batch-size", "3", "--sep", ";", "--encoding", "latin1"],
    )
    assert rc == 0
    assert called["batch_size"] == 3
    assert called["sep"] == ";"
    assert called["encoding"] == "latin1"
