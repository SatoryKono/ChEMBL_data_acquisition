"""Concurrency tests for :mod:`library.pipelines.testitem`."""

from __future__ import annotations

import threading
from pathlib import Path

import pandas as pd
import pytest

from types import SimpleNamespace

from library.config import ApiCfg, Config, PubChemCfg, RetryCfg
import library.pipelines.testitem as pipeline


def test_augment_pubchem_single_initialisation(monkeypatch) -> None:
    """Ensure concurrent augmentation initialises the PubChem session once."""

    threads = 5
    barrier = threading.Barrier(threads)
    init_lock = threading.Lock()
    init_calls = 0

    def fake_signature(api_cfg: ApiCfg, retry_cfg: RetryCfg) -> str:
        barrier.wait()
        return "signature"

    def fake_init_session(api_cfg: ApiCfg, retry_cfg: RetryCfg) -> None:
        nonlocal init_calls
        with init_lock:
            init_calls += 1

    def fake_add_pubchem_data(df: pd.DataFrame, *_args, **_kwargs) -> pd.DataFrame:
        return df

    monkeypatch.setattr(pipeline, "_PUBCHEM_SESSION_SIGNATURE", None)
    monkeypatch.setattr(pipeline, "_PUBCHEM_SESSION_LOCK", threading.Lock())
    monkeypatch.setattr(pipeline, "_pubchem_session_signature", fake_signature)
    monkeypatch.setattr(pipeline.pl, "init_session", fake_init_session)
    monkeypatch.setattr(
        pipeline, "_load_pubchem_cid_cache", lambda *_args, **_kwargs: {}
    )
    monkeypatch.setattr(pipeline, "add_pubchem_data", fake_add_pubchem_data)

    api_cfg = ApiCfg()
    retry_cfg = RetryCfg()
    pubchem_cfg = PubChemCfg()
    df = pd.DataFrame()
    errors: list[BaseException] = []

    def worker() -> None:
        try:
            pipeline.augment_pubchem(
                df,
                pubchem_cfg=pubchem_cfg,
                api_cfg=api_cfg,
                retry_cfg=retry_cfg,
                timeout=1.0,
                client=object(),
                fields=None,
                request_limit=10,
            )
        except BaseException as exc:  # pragma: no cover - defensive
            errors.append(exc)

    workers = [threading.Thread(target=worker) for _ in range(threads)]
    for thread in workers:
        thread.start()
    for thread in workers:
        thread.join()

    assert not errors
    assert init_calls == 1


def test_run_pipeline_streams_chunks(monkeypatch, tmp_path, cfg) -> None:
    """Ensure :func:`run_testitem_pipeline` processes chunks sequentially."""

    chunk_a = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]})
    chunk_b = pd.DataFrame({"molecule_chembl_id": ["CHEMBL3"]})
    captured_prepare: list[pd.DataFrame] = []
    captured_run: list[pd.DataFrame] = []
    captured_chunks: list[pd.DataFrame] = []

    def fake_read_ids(*_args, **_kwargs):
        return 0, pipeline.ReadInputIdsResult(
            ids_iter=iter(["CHEMBL1", "CHEMBL2", "CHEMBL3"]),
            sample_ids=("CHEMBL1",),
        )

    def fake_fetch(*_args, **_kwargs):
        return 0, iter([chunk_a, chunk_b]), ("CHEMBL1", "CHEMBL2", "CHEMBL3")

    def fake_prepare(df: pd.DataFrame, **_kwargs):
        captured_prepare.append(df.copy())
        lookup = pipeline.ParentLookupPreparedData(
            child_ids=pd.Series(dtype="string"),
            existing_parent_ids=pd.Series(dtype="string"),
            need_lookup=set(),
        )
        prep = pipeline.ParentEnrichmentPreparation(
            df=df,
            lookup_data=lookup,
            parent_catalog=None,
            parent_catalog_source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
            parent_stats=pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
                missing=0,
                unique=len(df),
                attached=len(df),
                uncovered=0,
            ),
        )
        return 0, prep

    def fake_run(prep: pipeline.ParentEnrichmentPreparation, **_kwargs):
        captured_run.append(prep.df.copy())
        stats = pipeline.ParentLookupStats(
            source=pipeline.PARENT_LOOKUP_SOURCE_LOOKUP,
            missing=0,
            unique=len(prep.df),
            attached=len(prep.df),
            uncovered=0,
        )
        return 0, pipeline.ParentEnrichmentResult(df=prep.df, parent_stats=stats)

    def fake_augment(df: pd.DataFrame, **_kwargs) -> pd.DataFrame:
        return df.assign(augmented=True)

    def fake_enrich(df: pd.DataFrame, **_kwargs):
        return 0, df

    def fake_finalize(chunks, **_kwargs):
        for chunk in chunks:
            captured_chunks.append(chunk.copy())
        return 0

    class DummyClient:
        def __init__(self, *_args, **_kwargs) -> None:  # pragma: no cover - simple stub
            self.client = SimpleNamespace()

        def __enter__(self) -> SimpleNamespace:
            return self.client

        def __exit__(self, *_args) -> bool:
            return False

    monkeypatch.setattr("library.pipelines.testitem.cli.read_input_ids", fake_read_ids)
    monkeypatch.setattr("library.pipelines.testitem.cli.fetch_testitems", fake_fetch)
    monkeypatch.setattr("library.pipelines.testitem.cli.prepare_parent_enrichment", fake_prepare)
    monkeypatch.setattr("library.pipelines.testitem.cli.run_parent_enrichment", fake_run)
    monkeypatch.setattr("library.pipelines.testitem.cli.augment_pubchem", fake_augment)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.apply_testitem_enrichment", fake_enrich
    )
    monkeypatch.setattr("library.pipelines.testitem.cli.finalize_output", fake_finalize)
    monkeypatch.setattr("library.pipelines.testitem.cli.ChemblClient", DummyClient)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.pc.init_session",
        lambda *_args, **_kwargs: None,
    )

    options = pipeline.TestitemPipelineOptions(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "output.csv",
    )

    result = pipeline.run_testitem_pipeline(cfg, options)

    assert result == 0
    assert [list(df["molecule_chembl_id"]) for df in captured_prepare] == [
        ["CHEMBL1", "CHEMBL2"],
        ["CHEMBL3"],
    ]
    assert [list(df["molecule_chembl_id"]) for df in captured_run] == [
        ["CHEMBL1", "CHEMBL2"],
        ["CHEMBL3"],
    ]
    assert [list(df["molecule_chembl_id"]) for df in captured_chunks] == [
        ["CHEMBL1", "CHEMBL2"],
        ["CHEMBL3"],
    ]


def test_run_pipeline_streams_missing(monkeypatch, tmp_path, cfg) -> None:
    """Missing identifiers are appended as a final placeholder chunk."""

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})
    captured_chunks: list[pd.DataFrame] = []

    def fake_read_ids(*_args, **_kwargs):
        return 0, pipeline.ReadInputIdsResult(
            ids_iter=iter(["CHEMBL1", "CHEMBL2"]),
            sample_ids=("CHEMBL1",),
        )

    def fake_fetch(*_args, **_kwargs):
        return 0, iter([chunk]), ("CHEMBL1", "CHEMBL2")

    def fake_prepare(df: pd.DataFrame, **_kwargs):
        lookup = pipeline.ParentLookupPreparedData(
            child_ids=pd.Series(dtype="string"),
            existing_parent_ids=pd.Series(dtype="string"),
            need_lookup=set(),
        )
        prep = pipeline.ParentEnrichmentPreparation(
            df=df,
            lookup_data=lookup,
            parent_catalog=None,
            parent_catalog_source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
            parent_stats=pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
                missing=0,
                unique=len(df),
                attached=len(df),
                uncovered=0,
            ),
        )
        return 0, prep

    def fake_run(prep: pipeline.ParentEnrichmentPreparation, **_kwargs):
        stats = pipeline.ParentLookupStats(
            source=pipeline.PARENT_LOOKUP_SOURCE_LOOKUP,
            missing=0,
            unique=len(prep.df),
            attached=len(prep.df),
            uncovered=0,
        )
        return 0, pipeline.ParentEnrichmentResult(df=prep.df, parent_stats=stats)

    def fake_finalize(chunks, **_kwargs):
        for item in chunks:
            captured_chunks.append(item.copy())
        return 0

    class DummyClient:
        def __init__(self, *_args, **_kwargs) -> None:
            self.client = SimpleNamespace()

        def __enter__(self) -> SimpleNamespace:
            return self.client

        def __exit__(self, *_args) -> bool:
            return False

    monkeypatch.setattr("library.pipelines.testitem.cli.read_input_ids", fake_read_ids)
    monkeypatch.setattr("library.pipelines.testitem.cli.fetch_testitems", fake_fetch)
    monkeypatch.setattr("library.pipelines.testitem.cli.prepare_parent_enrichment", fake_prepare)
    monkeypatch.setattr("library.pipelines.testitem.cli.run_parent_enrichment", fake_run)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.augment_pubchem", lambda df, **_: df
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.apply_testitem_enrichment",
        lambda df, **_: (0, df),
    )
    monkeypatch.setattr("library.pipelines.testitem.cli.finalize_output", fake_finalize)
    monkeypatch.setattr("library.pipelines.testitem.cli.ChemblClient", DummyClient)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.pc.init_session",
        lambda *_args, **_kwargs: None,
    )

    options = pipeline.TestitemPipelineOptions(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "output.csv",
    )

    result = pipeline.run_testitem_pipeline(cfg, options)

    assert result == 0
    assert [list(df["molecule_chembl_id"]) for df in captured_chunks] == [
        ["CHEMBL1"],
        ["CHEMBL2"],
    ]


def test_run_pipeline_parentless_skips_full_sync(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, cfg: Config
) -> None:
    """Legacy wrapper also honours parentless full sync flag."""

    cfg.molecule_catalog.enable_full_sync_for_parentless = False
    cfg.molecule_catalog.hierarchy_lookup_path = None
    cfg.molecule_catalog.cache_path = tmp_path / "catalog.json"
    cfg.molecule_catalog.sqlite_path = tmp_path / "catalog.sqlite"
    cfg.pubchem.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            cfg.molecule_catalog.parent_field: [pd.NA, pd.NA],
        }
    )

    def fake_read_ids(*_args, **_kwargs):
        return 0, pipeline.ReadInputIdsResult(
            ids_iter=iter(["CHEMBL1", "CHEMBL2"]),
            sample_ids=("CHEMBL1",),
        )

    def fake_fetch(*_args, **_kwargs):
        return 0, iter([chunk]), ("CHEMBL1", "CHEMBL2")

    captured_stats: dict[str, pipeline.ParentLookupStats] = {}

    def fake_finalize(chunks, **kwargs):
        for _chunk in chunks:
            pass
        captured_stats["parent"] = kwargs["parent_stats_supplier"]()
        return 0

    def fail_full_sync(*_args, **_kwargs):
        raise AssertionError("full sync should be disabled for parentless identifiers")

    class DummyClient:
        def __init__(self, *_args, **_kwargs) -> None:
            self.client = SimpleNamespace()

        def __enter__(self) -> SimpleNamespace:
            return self.client

        def __exit__(self, *_args) -> bool:
            return False

    monkeypatch.setattr("library.testitem_pipeline.cli.read_input_ids", fake_read_ids)
    monkeypatch.setattr("library.testitem_pipeline.cli.fetch_testitems", fake_fetch)
    monkeypatch.setattr("library.testitem_pipeline.cli.finalize_output", fake_finalize)
    monkeypatch.setattr("library.testitem_pipeline.cli.augment_pubchem", lambda df, **_: df)
    monkeypatch.setattr(
        "library.testitem_pipeline.cli.apply_testitem_enrichment", lambda df, **_: (0, df)
    )
    monkeypatch.setattr("library.testitem_pipeline.cli.ChemblClient", DummyClient)
    monkeypatch.setattr(
        "library.pipelines.testitem.catalog.molecule_catalog.fetch_parent_catalog_for",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.catalog.load_parent_catalog",
        fail_full_sync,
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.catalog.query_parent_catalog", lambda *_args, **_kwargs: {}
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.catalog.update_parent_catalog_cache",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.catalog.write_parent_catalog_cache",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "library.testitem_pipeline.catalog.molecule_catalog.fetch_parent_catalog_for",
        lambda *_args, **_kwargs: {},
    )
    monkeypatch.setattr(
        "library.testitem_pipeline.catalog.load_parent_catalog",
        fail_full_sync,
    )
    monkeypatch.setattr(
        "library.testitem_pipeline.catalog.query_parent_catalog", lambda *_args, **_kwargs: {}
    )
    monkeypatch.setattr(
        "library.testitem_pipeline.catalog.update_parent_catalog_cache",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        "library.testitem_pipeline.catalog.write_parent_catalog_cache",
        lambda *_args, **_kwargs: None,
    )

    options = pipeline.TestitemPipelineOptions(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "output.csv",
    )

    result = pipeline.run_testitem_pipeline(cfg, options)

    assert result == 0
    parent_stats = captured_stats.get("parent")
    assert parent_stats is not None
    assert parent_stats.full_sync_duration_seconds == 0.0
    assert parent_stats.remaining_before_full_sync == 0

def test_run_pipeline_skips_pubchem_when_disabled(monkeypatch, tmp_path, cfg) -> None:
    """Ensure PubChem augmentation is bypassed when disabled in configuration."""

    cfg.pubchem.enable = False
    placeholder = "chembl-da/0.1 (mailto:contact@example.org)"
    cfg.api.user_agent = placeholder
    cfg.pubchem.user_agent = placeholder

    chunk = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    def fake_read_ids(*_args, **_kwargs):
        return 0, pipeline.ReadInputIdsResult(
            ids_iter=iter(["CHEMBL1"]),
            sample_ids=("CHEMBL1",),
        )

    def fake_fetch(*_args, **_kwargs):
        return 0, iter([chunk]), ("CHEMBL1",)

    def fake_prepare(df: pd.DataFrame, **_kwargs):
        lookup = pipeline.ParentLookupPreparedData(
            child_ids=pd.Series(dtype="string"),
            existing_parent_ids=pd.Series(dtype="string"),
            need_lookup=set(),
        )
        prep = pipeline.ParentEnrichmentPreparation(
            df=df,
            lookup_data=lookup,
            parent_catalog=None,
            parent_catalog_source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
            parent_stats=pipeline.ParentLookupStats(
                source=pipeline.PARENT_LOOKUP_SOURCE_CACHE,
                missing=0,
                unique=len(df),
                attached=len(df),
                uncovered=0,
            ),
        )
        return 0, prep

    def fake_run(prep: pipeline.ParentEnrichmentPreparation, **_kwargs):
        stats = pipeline.ParentLookupStats(
            source=pipeline.PARENT_LOOKUP_SOURCE_LOOKUP,
            missing=0,
            unique=len(prep.df),
            attached=len(prep.df),
            uncovered=0,
        )
        return 0, pipeline.ParentEnrichmentResult(df=prep.df, parent_stats=stats)

    def fake_enrich(df: pd.DataFrame, **_kwargs):
        return 0, df

    def fake_finalize(*_args, **_kwargs):
        return 0

    def fail_prepare_pubchem(*_args, **_kwargs):
        raise AssertionError("PubChem configuration should not be prepared when disabled")

    def fail_init_session(*_args, **_kwargs):
        raise AssertionError("PubChem session should not be initialised when disabled")

    def fail_augment(*_args, **_kwargs):
        raise AssertionError("augment_pubchem should not run when disabled")

    class DummyClient:
        def __init__(self, *_args, **_kwargs) -> None:  # pragma: no cover - simple stub
            self.client = SimpleNamespace()

        def __enter__(self) -> SimpleNamespace:
            return self.client

        def __exit__(self, *_args) -> bool:
            return False

    monkeypatch.setattr("library.pipelines.testitem.cli.read_input_ids", fake_read_ids)
    monkeypatch.setattr("library.pipelines.testitem.cli.fetch_testitems", fake_fetch)
    monkeypatch.setattr("library.pipelines.testitem.cli.prepare_parent_enrichment", fake_prepare)
    monkeypatch.setattr("library.pipelines.testitem.cli.run_parent_enrichment", fake_run)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.apply_testitem_enrichment", fake_enrich
    )
    monkeypatch.setattr("library.pipelines.testitem.cli.finalize_output", fake_finalize)
    monkeypatch.setattr("library.pipelines.testitem.cli.ChemblClient", DummyClient)
    monkeypatch.setattr(
        "library.pipelines.testitem.cli._prepare_pubchem_api_cfg",
        fail_prepare_pubchem,
    )
    monkeypatch.setattr(
        "library.pipelines.testitem.cli.pc.init_session", fail_init_session
    )
    monkeypatch.setattr("library.pipelines.testitem.cli.augment_pubchem", fail_augment)

    options = pipeline.TestitemPipelineOptions(
        input_csv=tmp_path / "input.csv",
        output_csv=tmp_path / "output.csv",
    )

    result = pipeline.run_testitem_pipeline(cfg, options)

    assert result == 0

