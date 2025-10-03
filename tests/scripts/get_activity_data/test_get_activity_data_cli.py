"""Smoke tests for :mod:`scripts.get_activity_data`."""

from importlib import import_module
from pathlib import Path
import time

import pandas as pd
import pytest
import yaml

from library import cli_utils


def test_get_activity_data_cli_exposes_main() -> None:
    module = import_module("scripts.get_activity_data")
    assert hasattr(module, "main")


def test_get_activity_data_cli_print_config() -> None:
    module = import_module("scripts.get_activity_data")
    # ``--print-config`` ensures the pipeline stops before network or IO work.
    exit_code = module.main(["--print-config"])
    assert exit_code == 0


def test_get_activity_data_cli_workers_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = import_module("scripts.get_activity_data")

    output_dir = tmp_path / "output"
    cache_dir = tmp_path / "cache"
    config_path = tmp_path / "config.yaml"
    input_path = tmp_path / "activity.csv"
    output_path = tmp_path / "activities.csv"

    config_data = {
        "sources": {
            "chembl": {
                "pipelines": {
                    "activity": {
                        "workers": 2,
                        "batch_size": 2,
                    }
                }
            }
        },
        "local": {
            "io": {
                "output_dir": str(output_dir),
                "cache_dir": str(cache_dir),
                "exist_ok": True,
            }
        },
    }

    config_path.write_text(yaml.safe_dump(config_data))
    input_path.write_text("activity_id\nA1\nA2\nA3\nA4\n")

    monkeypatch.setattr(cli_utils, "ensure_dirs", lambda cfg: None)
    monkeypatch.setattr(
        module.io,
        "read_ids",
        lambda *args, **kwargs: iter(["A1", "A2", "A3", "A4"]),
    )

    captured = {"chunks": [], "workers": None, "calls": []}

    frames = {
        ("A1", "A2"): pd.DataFrame({"activity_id": ["chunk0_A1", "chunk0_A2"]}),
        ("A3", "A4"): pd.DataFrame({"activity_id": ["chunk1_A3", "chunk1_A4"]}),
    }
    delays = {("A1", "A2"): 0.05, ("A3", "A4"): 0.0}

    def fake_get_activities(ids, *, cfg, client, chunk_size, timeout, **kwargs):
        key = tuple(ids)
        captured["calls"].append((key, chunk_size, timeout))
        delay = delays.get(key, 0.0)
        if delay:
            time.sleep(delay)
        return frames[key]

    class DummyClient:
        def __init__(self, *args, **kwargs) -> None:  # pragma: no cover - trivial
            pass

        def __enter__(self):  # pragma: no cover - trivial
            return self

        def __exit__(self, *args):  # pragma: no cover - trivial
            return False

    monkeypatch.setattr(module.cl, "get_activities", fake_get_activities)
    monkeypatch.setattr(module, "ChemblClient", lambda *args, **kwargs: DummyClient())
    monkeypatch.setattr(module, "get_global_limiter", lambda *args, **kwargs: None)

    def _normalize_chunks(source):
        if source is None:
            return []
        if isinstance(source, pd.DataFrame):
            return [source]
        return list(source)

    def fake_run_pipeline(*, fetcher, cfg, **kwargs):
        captured["workers"] = cfg.activity.workers
        iterator = fetcher()
        collected: list[pd.DataFrame] = []
        while True:
            try:
                chunk = next(iterator)
            except StopIteration as stop:
                collected.extend(_normalize_chunks(stop.value))
                break
            else:
                collected.append(chunk)
        for frame in collected:
            captured["chunks"].append(list(frame["activity_id"]))
        return 0

    monkeypatch.setattr(module, "run_pipeline", fake_run_pipeline)

    exit_code = module.main(
        [
            "--config",
            str(config_path),
            "--input",
            str(input_path),
            "--final-out",
            str(output_path),
        ]
    )

    assert exit_code == 0
    assert captured["workers"] == 2
    assert captured["chunks"] == [
        ["chunk0_A1", "chunk0_A2"],
        ["chunk1_A3", "chunk1_A4"],
    ]
    assert captured["calls"] == [
        (("A1", "A2"), 2, 30.0),
        (("A3", "A4"), 2, 30.0),
    ]
