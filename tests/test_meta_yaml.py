from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

from library.config import Config
from library.io import write_csv
from library.sidecar import SidecarErrors
from schemas import CsvMetaSchema


def test_write_csv_creates_meta(tmp_path: Path) -> None:
    df = pd.DataFrame({"id": [1, 2], "name": ["a", "b"]})
    cfg = Config()
    out_path = tmp_path / "out.csv"
    write_csv(df, out_path, cfg=cfg, key_cols=["id"])
    meta_path = Path(f"{out_path}.meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf8"))
    CsvMetaSchema.model_validate(meta)
    assert meta["columns"] == ["id", "name"]
    assert meta["dtypes"] == {"id": "int64", "name": "object"}
    assert meta["git_sha"]
    datetime.fromisoformat(meta["generated_at"])  # raises if invalid


def test_sidecar_writes_meta(tmp_path: Path) -> None:
    errors = SidecarErrors()
    errors.add_error({"col_a": "foo", "col_b": 1})
    out = tmp_path / "errors.csv"
    errors.save(out, cfg=Config())
    meta_path = Path(f"{out}.meta.yaml")
    assert meta_path.exists()
    meta = yaml.safe_load(meta_path.read_text(encoding="utf8"))
    CsvMetaSchema.model_validate(meta)
    assert meta["columns"] == ["col_a", "col_b"]
    assert meta["dtypes"] == {"col_a": "string", "col_b": "string"}
    assert meta["git_sha"]
    datetime.fromisoformat(meta["generated_at"])  # raises if invalid
