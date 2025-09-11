from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml

from library import io
from library.config import Config, IoCfg


def test_read_csv_validates_columns(tmp_path: Path) -> None:
    path = tmp_path / "data.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["a", "b"])
        writer.writerow(["1", "2"])
    df = io.read_csv(path, cfg=IoCfg(), required_columns=["a"])
    assert list(df.columns) == ["a", "b"]


def test_read_csv_missing_column(tmp_path: Path) -> None:
    path = tmp_path / "data.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["a"])
        writer.writerow(["1"])
    with pytest.raises(ValueError):
        io.read_csv(path, cfg=IoCfg(), required_columns=["a", "b"])


def test_write_csv_creates_single_metadata_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    df = pd.DataFrame({"a": [1]})
    path = tmp_path / "out.csv"
    cfg = Config()

    calls: list[Path] = []

    def spy(p: Path, cfg: Config) -> None:
        calls.append(p)
        # Create a dummy metadata file to emulate the original behaviour
        Path(str(p) + ".meta.yaml").touch()

    monkeypatch.setattr(io, "_write_meta", spy)
    io.write_csv(df, path, cfg=cfg)

    meta_path = Path(str(path) + ".meta.yaml")
    assert len(calls) == 1
    assert path.exists()
    assert meta_path.exists()


def test_write_meta_serialises_paths(tmp_path: Path) -> None:
    df = pd.DataFrame({"a": [1]})
    path = tmp_path / "out.csv"
    cfg = Config()
    io.write_csv(df, path, cfg=cfg)

    meta_path = Path(f"{path}.meta.yaml")
    meta = yaml.safe_load(meta_path.read_text())
    assert isinstance(meta["config"]["io"]["output_dir"], str)
    assert meta["config"]["io"]["output_dir"] == str(cfg.io.output_dir)


def test_write_csv_deterministic_hash(tmp_path: Path) -> None:
    df = pd.DataFrame({"b": [3, 1], "a": [4.0, 2.0]})
    path = tmp_path / "out.csv"
    cfg = Config()
    io.write_csv(df, path, cfg=cfg)
    first_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    shuffled = df.sample(frac=1).reset_index(drop=True)[["b", "a"]]
    io.write_csv(shuffled, path, cfg=cfg)
    second_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert first_hash == second_hash
