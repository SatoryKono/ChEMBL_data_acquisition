from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from io import StringIO
from pathlib import Path

import pandas as pd
import pandera as pa
import pytest
import yaml

import library.git_utils as git_utils
from library import io
from library.config import Config, IoCfg
from library.utils.logging_setup import LoggerConfig, configure_logger


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


def test_read_ids_missing_column_lists_available(tmp_path: Path) -> None:
    path = tmp_path / "data.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["a", "b"])
        writer.writerow(["1", "2"])
    with pytest.raises(ValueError) as exc:
        list(io.read_ids(path, column="c", cfg=IoCfg()))
    assert (
        str(exc.value)
        == f"column 'c' not found in {path}; available columns: ['a', 'b']"
    )


def test_read_csv_types_and_na_values() -> None:
    path = Path("tests/data/io_types.csv")
    schema = pa.DataFrameSchema(
        {
            "flag": pa.Column("boolean", nullable=True),
            "date": pa.Column("datetime64[ns]", nullable=True),
            "count": pa.Column("Int64", nullable=True),
        },
    )
    df = io.read_csv(
        path,
        cfg=IoCfg(),
        dtype={"flag": "boolean", "count": "Int64"},
        na_values=["---", "missing", ""],
        parse_dates=["date"],
        schema=schema,
    )
    assert df["flag"].dtype == "boolean"
    assert df["date"].dtype == "datetime64[ns]"
    assert df["count"].dtype == "Int64"
    assert pd.isna(df.loc[2, "flag"])  # custom NA token
    assert pd.isna(df.loc[1, "date"])  # empty string -> NA
    assert pd.isna(df.loc[2, "count"])  # custom NA token


def test_read_csv_with_schema(tmp_path: Path) -> None:
    """``read_csv`` validates against a provided schema."""

    path = tmp_path / "data.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["a", "b"])
        writer.writerow(["1", "x"])

    good_schema = pa.DataFrameSchema({"a": pa.Column(int), "b": pa.Column(str)})
    df = io.read_csv(path, cfg=IoCfg(), schema=good_schema)
    assert list(df.columns) == ["a", "b"]

    bad_schema = pa.DataFrameSchema({"a": pa.Column(int), "c": pa.Column(str)})
    with pytest.raises(pa.errors.SchemaError):
        io.read_csv(path, cfg=IoCfg(), schema=bad_schema)


def test_write_csv_creates_metadata_file(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [1]})
    path = tmp_path / "out.csv"
    io.write_csv(df, path, cfg=cfg)
    meta_files = list(tmp_path.glob("*.meta.yaml"))
    assert path.exists()
    assert len(meta_files) == 1
    assert meta_files[0].exists()


def test_write_meta_serialises_paths(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [1]})
    path = tmp_path / "out.csv"
    io.write_csv(df, path, cfg=cfg)

    meta_path = Path(f"{path}.meta.yaml")
    meta = yaml.safe_load(meta_path.read_text())
    assert isinstance(meta["config"]["local"]["io"]["output_dir"], str)
    assert meta["config"]["local"]["io"]["output_dir"] == str(cfg.io.output_dir)


def test_write_csv_deterministic_hash(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"b": [3, 1], "a": [4.0, 2.0]})
    path = tmp_path / "out.csv"
    io.write_csv(df, path, cfg=cfg)
    first_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    shuffled = df.sample(frac=1).reset_index(drop=True)[["b", "a"]]
    io.write_csv(shuffled, path, cfg=cfg)
    second_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert first_hash == second_hash


def test_write_csv_with_key_cols(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [2, 1], "b": ["x", "y"]})
    path = tmp_path / "out.csv"
    io.write_csv(df, path, cfg=cfg, key_cols=["a"])
    first_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    shuffled = df.sample(frac=1).reset_index(drop=True)
    io.write_csv(shuffled, path, cfg=cfg, key_cols=["a"])
    second_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    assert first_hash == second_hash


def test_read_ids_custom_na_marker(tmp_path: Path) -> None:
    """``read_ids`` honours custom NA markers."""

    path = tmp_path / "ids.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["id"])
        writer.writerow(["1"])
        writer.writerow(["NA"])
        writer.writerow(["2"])

    ids = list(io.read_ids(path, column="id", cfg=IoCfg(), na_markers=["NA"]))
    assert ids == ["1", "2"]


def test_read_ids_falls_back_to_alternative_encoding(tmp_path: Path) -> None:
    """``read_ids`` retries with configured fallback encodings."""

    path = tmp_path / "ids.csv"
    path.write_bytes("id\nCHEMBL±1\n".encode("windows-1251"))

    cfg = IoCfg(csv_encoding="utf-8-sig")
    ids = list(io.read_ids(path, column="id", cfg=cfg))
    assert ids == ["CHEMBL±1"]


def test_read_ids_uses_locale_encoding_when_config_lacks_fallback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Locale preferred encoding is appended when custom fallbacks are absent."""

    path = tmp_path / "ids.csv"
    path.write_bytes("id\nЖ\n".encode("windows-1251"))

    cfg = IoCfg(csv_encoding="utf-8", csv_fallback_encodings=())

    monkeypatch.setattr(io.locale, "getpreferredencoding", lambda _=False: "windows-1251")

    ids = list(io.read_ids(path, column="id", cfg=cfg))
    assert ids == ["Ж"]


def test_read_ids_raises_when_all_encodings_fail(tmp_path: Path) -> None:
    """``read_ids`` surfaces decoding errors after exhausting fallbacks."""

    path = tmp_path / "ids.csv"
    path.write_bytes("id\n1\n".encode("utf-16"))

    cfg = IoCfg(csv_encoding="utf-8", csv_fallback_encodings=())
    with pytest.raises(ValueError) as exc:
        list(io.read_ids(path, column="id", cfg=cfg))
    assert "failed to decode CSV" in str(exc.value)


def test_write_csv_missing_key_column(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame({"a": [1], "b": [2]})
    path = tmp_path / "out.csv"
    with pytest.raises(ValueError, match="Missing key columns: \['c'\]"):
        io.write_csv(df, path, cfg=cfg, key_cols=["c"])


def test_write_csv_normalises_types(tmp_path: Path, cfg: Config) -> None:
    df = pd.DataFrame(
        {
            "b": [True, False],
            "d": [pd.Timestamp("2020-01-02"), pd.Timestamp("2020-01-01")],
        }
    )
    path = tmp_path / "out.csv"
    io.write_csv(df, path, cfg=cfg, key_cols=["d"])
    text = path.read_text(encoding="utf-8-sig")
    assert text == ("b,d\nfalse,2020-01-01\ntrue,2020-01-02\n")


def test_write_csv_chunksize(tmp_path: Path, cfg: Config) -> None:
    """Chunked writes yield the same output as unchunked writes."""

    rows = 1500
    df = pd.DataFrame(
        {
            "a": range(rows),
            "b": [i % 2 == 0 for i in range(rows)],
            "d": pd.date_range("2020-01-01", periods=rows, freq="D"),
            "f": [i / 3 for i in range(rows)],
        }
    )
    path1 = tmp_path / "plain.csv"
    path2 = tmp_path / "chunked.csv"
    io.write_csv(df, path1, cfg=cfg)
    io.write_csv(df, path2, cfg=cfg, chunksize=500)
    assert path1.read_bytes() == path2.read_bytes()


def test_git_sha_timeout_returns_unknown(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """_git_sha returns ``"unknown"`` when the git command exceeds timeout."""

    stream = StringIO()
    monkeypatch.setattr(
        git_utils,
        "logger",
        configure_logger(LoggerConfig(level="WARNING", stream=stream)),
    )

    def fail(*args: object, **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise subprocess.CalledProcessError(
            returncode=1,
            cmd=["git", "rev-parse", "HEAD"],
            stderr="fatal: simulated error",
        )

    monkeypatch.setattr(git_utils.subprocess, "run", fail)
    monkeypatch.setattr(git_utils, "_read_head_sha", lambda *_: None)
    git_utils._git_sha.cache_clear()

    assert git_utils._git_sha() == "UNKNOWN"
    record = json.loads(stream.getvalue().splitlines()[0])
    assert record["event"] == "git_sha_unavailable"
    assert "error" in record
    assert record["error_returncode"] == 1
    assert record["error_returncode_raw"] == 1
    assert record["error_cmd"] == ["git", "rev-parse", "HEAD"]
    assert record["error_stderr"] == "fatal: simulated error"
