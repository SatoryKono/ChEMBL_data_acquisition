"""Concurrency tests for :mod:`library.testitem_pipeline`."""

from __future__ import annotations

import threading

import pandas as pd

from library.config import ApiCfg, PubChemCfg, RetryCfg
from library import testitem_pipeline as pipeline


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
    monkeypatch.setattr(pipeline, "_load_pubchem_cid_cache", lambda *_args, **_kwargs: {})
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
