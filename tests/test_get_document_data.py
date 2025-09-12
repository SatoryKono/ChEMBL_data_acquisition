from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from library import chembl_library as cl
from library import io
from scripts import get_document_data as gdd


def test_cli_uses_custom_column(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure CLI respects the ``--column`` option for ChEMBL documents."""
    input_csv = tmp_path / "docs.csv"
    input_csv.write_text("document_chembl_id\nCHEMBL1\n")

    captured: dict[str, str] = {}
    orig_read_ids = io.read_ids

    def fake_read_ids(
        path: str | Path,
        *,
        column: str,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
    ) -> Iterable[str]:
        captured["column"] = column
        return orig_read_ids(path, column=column, cfg=cfg, sep=sep, encoding=encoding)

    monkeypatch.setattr(io, "read_ids", fake_read_ids)
    monkeypatch.setattr(gdd, "ChemblClient", lambda *_, **__: object())

    def fake_get_documents(
        ids: Iterable[str],
        cfg: Any,
        client: Any,
        chunk_size: int,
        timeout: float,
    ) -> pd.DataFrame:
        data = list(ids)
        return pd.DataFrame({"document_chembl_id": data, "title": ["t"] * len(data)})

    monkeypatch.setattr(cl, "get_documents", fake_get_documents)
    monkeypatch.setattr(gdd, "normalize_documents", lambda df: df)

    def fake_write_csv(
        df: pd.DataFrame,
        path: Path,
        *,
        cfg: Any,
        sep: str | None = None,
        encoding: str | None = None,
        key_cols: list[str] | None = None,
        chunksize: int | None = None,
    ) -> Path:
        return path

    monkeypatch.setattr(io, "write_csv", fake_write_csv)
    monkeypatch.setattr(gdd, "file_sha256", lambda p: "deadbeef")
    monkeypatch.setattr(gdd, "write_meta_yaml", lambda **__: None)
    monkeypatch.setattr(gdd, "analyze_table_quality", lambda df, table_name: None)
    monkeypatch.setattr(gdd, "ensure_dirs", lambda cfg: None)

    rc = gdd.main(
        [
            "chembl",
            "--column",
            "document_chembl_id",
            "--input",
            str(input_csv),
            "--output",
            str(tmp_path / "out.csv"),
        ]
    )
    assert rc == 0
    assert captured["column"] == "document_chembl_id"
