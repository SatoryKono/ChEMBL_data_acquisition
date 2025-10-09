"""Profile the get_activity_data pipeline with deterministic stubs."""

from __future__ import annotations

# ruff: noqa: E402  # requires sys.path mutation before local imports

import argparse
import math
import sys
from collections.abc import Callable
from contextlib import ExitStack
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from library.config import Config


RESOURCE_DIR = ROOT / "tests" / "resources" / "activity_pipeline"


class StageTimer:
    def __init__(self) -> None:
        self.durations: dict[str, float] = {
            "load": 0.0,
            "fetch": 0.0,
            "process": 0.0,
            "postprocess": 0.0,
            "save": 0.0,
        }

    def wrap(self, stage: str, func: Callable[..., object]) -> Callable[..., object]:
        def wrapper(*args, **kwargs):
            start = perf_counter()
            try:
                return func(*args, **kwargs)
            finally:
                self.durations[stage] = self.durations.get(stage, 0.0) + perf_counter() - start

        return wrapper


def profile(limit: int | None = None, *, rows: int = 300) -> dict[str, float]:
    from scripts import get_activity_data

    base_chunk = pd.read_csv(RESOURCE_DIR / "chunk_happy.csv")
    repeats = max(1, math.ceil(rows / len(base_chunk)))
    chunk_df = pd.concat([base_chunk] * repeats, ignore_index=True)
    chunk_df = chunk_df.iloc[:rows].copy()
    chunk_df["activity_id"] = [f"ACT{i+1}" for i in range(len(chunk_df))]
    chunk_df["molecule_chembl_id"] = [f"CHEMBL{i+1}" for i in range(len(chunk_df))]
    chunk_df["assay_chembl_id"] = [f"ASSAY{i+1}" for i in range(len(chunk_df))]

    timers = StageTimer()
    stack = ExitStack()

    mock = __import__("unittest.mock", fromlist=["patch"])
    stack.enter_context(mock.patch("library.resources.dictionaries._parse_manifest", return_value={}))
    stack.enter_context(
        mock.patch(
            "library.resources.dictionaries.get_resource",
            side_effect=lambda name, base_dir=None: SimpleNamespace(
                name=name,
                relative_path=Path(name),
                path=ROOT / "tests" / "resources" / "molecule_catalog.csv",
                version="test",
                sha256="dummy",
                generator=Path("generator.py"),
            ),
        )
    )
    stack.enter_context(
        mock.patch(
            "library.resources.dictionaries.get_resource_path",
            side_effect=lambda name, base_dir=None: ROOT
            / "tests"
            / "resources"
            / "molecule_catalog.csv",
        )
    )

    cfg = Config()
    cfg.activity.limit = limit if limit is not None else rows
    cfg.activity.dry_run = False
    cfg.activity.batch_size = 3
    cfg.activity.workers = 1
    cfg.system.doc_quality.enable = False
    cfg.activity_enrichment.action_type.enabled = False
    cfg.activity_enrichment.activity_properties.enabled = False

    testitem_frame = pd.DataFrame(
        {
            "molecule_chembl_id": [f"CHEMBL{i+1}" for i in range(len(chunk_df))],
            "pref_name": [f"MOLECULE_{i+1}" for i in range(len(chunk_df))],
        }
    )

    def fake_read_ids(path: Path, *, column: str, cfg):  # type: ignore[override]
        del cfg
        start = perf_counter()
        try:
            frame = pd.read_csv(path, dtype="string")
            return frame[column].astype("string").tolist()
        finally:
            timers.durations["load"] += perf_counter() - start

    def fake_get_activities(chunk_ids, **_kwargs):
        start = perf_counter()
        try:
            identifiers = [str(item) for item in chunk_ids]
            mask = chunk_df["activity_id"].astype("string").isin(identifiers)
            return chunk_df.loc[mask].reset_index(drop=True)
        finally:
            timers.durations["fetch"] += perf_counter() - start

    def fake_get_testitem(chunk_ids, **_kwargs):
        identifiers = [str(item) for item in chunk_ids]
        mask = testitem_frame["molecule_chembl_id"].astype("string").isin(identifiers)
        return testitem_frame.loc[mask].reset_index(drop=True)

    class DummyClient:
        def __init__(self, *_, **__):
            pass

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def close(self):
            return None

    def instrument(stage: str, qualname: str) -> None:
        module_name, _, attr = qualname.rpartition(".")
        module = __import__(module_name, fromlist=[attr])
        original = getattr(module, attr)
        setattr(module, attr, timers.wrap(stage, original))
        stack.callback(lambda: setattr(module, attr, original))

    stack.enter_context(
        mock.patch("scripts.get_activity_data.io.read_ids", side_effect=fake_read_ids)
    )
    stack.enter_context(
        mock.patch("scripts.get_activity_data.cl.get_activities", side_effect=fake_get_activities)
    )
    stack.enter_context(
        mock.patch("scripts.get_activity_data.cl.get_testitem", side_effect=fake_get_testitem)
    )
    try:
        stack.enter_context(mock.patch("library.orchestration.context.ChemblClient", DummyClient))
    except AttributeError:
        stack.enter_context(mock.patch("library.clients.ChemblClient", DummyClient))
    assay_lookup = {f"ASSAY{i+1}": f"SRC-ASSAY{i+1}" for i in range(len(chunk_df))}
    stack.enter_context(
        mock.patch(
            "scripts.get_activity_data._load_assay_src_lookup",
            return_value=assay_lookup,
        )
    )
    stack.enter_context(
        mock.patch(
            "scripts.get_activity_data.process_activity_extended",
            timers.wrap("postprocess", lambda *args, **kwargs: None),
        )
    )

    instrument("process", "scripts.get_activity_data.compute_activity_bounds")
    instrument("process", "scripts.get_activity_data.apply_activity_annotations")
    instrument("process", "scripts.get_activity_data.normalize_activities")
    instrument("process", "scripts.get_activity_data.add_pipeline_metadata")

    stack.enter_context(
        mock.patch(
            "scripts.get_activity_data.write_csv_chunks_deterministic",
            new=timers.wrap("save", get_activity_data.write_csv_chunks_deterministic),
        )
    )

    input_path = Path("input.csv")
    with input_path.open("w", encoding="utf-8") as handle:
        handle.write("activity_id\n")
        for idx in range(len(chunk_df)):
            handle.write(f"ACT{idx+1}\n")

    args = argparse.Namespace(
        input_csv=input_path,
        output_csv=Path("activities.csv"),
        final_out=Path("activities.csv"),
        offset=0,
        workers=None,
        skip_existing=False,
        force=False,
        dry_run=False,
        invocation=None,
    )

    try:
        get_activity_data.run(cfg, args)
    finally:
        stack.close()
        try:
            input_path.unlink()
        except FileNotFoundError:
            pass
        try:
            Path("activities.csv").unlink()
        except FileNotFoundError:
            pass

    return timers.durations


if __name__ == "__main__":
    import json

    durations = profile()
    print(json.dumps(durations, indent=2))
