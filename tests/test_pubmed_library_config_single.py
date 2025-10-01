from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library import pubmed_library as pl
from library.config import RetryCfg


class DummyLimiter:
    def acquire(self) -> None:  # pragma: no cover - simple stub
        return None


def _stub_network(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_fetch_pubmed_batch(
        session,
        pmids,
        delay,
        cfg=None,
        *,
        retry_cfg: RetryCfg | None = None,
        client=None,
    ):
        return [
            {
                "PubMed.PMID": pid,
                "PubMed.DOI": "10.1000/example",
                "PubMed.ArticleTitle": "Example",
            }
            for pid in pmids
        ]

    def fake_fetch_semantic_scholar_batch(session, pmids, delay, cfg=None):
        return [
            {
                "scholar.PMID": pid,
                "scholar.DOI": "10.1000/example",
                "scholar.Error": "",
            }
            for pid in pmids
        ]

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        pl, "fetch_semantic_scholar_batch", fake_fetch_semantic_scholar_batch
    )
    monkeypatch.setattr(pl, "fetch_openalex", lambda *a, **k: {"OpenAlex.Error": ""})
    monkeypatch.setattr(pl, "fetch_crossref", lambda *a, **k: {"crossref.Error": ""})
    monkeypatch.setattr(pl, "get_limiter", lambda *a, **k: DummyLimiter())


def test_pubmed_library_single_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure :func:`pubmed_library.main` creates ``Config`` only once."""
    _stub_network(monkeypatch)
    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "out.csv"

    calls: dict[str, int] = {"cfg": 0}
    OriginalConfig = pl.Config

    class CountingConfig(OriginalConfig):
        def __init__(self, *args, **kwargs):
            calls["cfg"] += 1
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(pl, "Config", CountingConfig)

    exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
        ]
    )
    assert exit_code == 0
    assert output_csv.exists()
    df = pd.read_csv(output_csv)
    assert not df.empty
    assert calls["cfg"] == 1


def test_pubmed_library_rate_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PubMed-specific rate overrides should control limiter and delays."""

    calls: dict[str, list[tuple[str, float, int | None]]] = {"limiters": []}
    delays: dict[str, float] = {}

    def fake_fetch_pubmed_batch(
        session,
        pmids,
        delay,
        cfg=None,
        *,
        retry_cfg: RetryCfg | None = None,
        client=None,
    ):
        delays["pubmed"] = delay
        return [
            {
                "PubMed.PMID": pid,
                "PubMed.DOI": "10.1000/example",  # pragma: no cover - deterministic field
                "PubMed.ArticleTitle": "Example",
            }
            for pid in pmids
        ]

    def fake_fetch_semantic_scholar_batch(session, pmids, delay, cfg=None):
        delays["semantic_scholar"] = delay
        return [
            {
                "scholar.PMID": pid,
                "scholar.DOI": "10.1000/example",  # pragma: no cover - deterministic field
                "scholar.Error": "",
            }
            for pid in pmids
        ]

    monkeypatch.setattr(pl, "fetch_pubmed_batch", fake_fetch_pubmed_batch)
    monkeypatch.setattr(
        pl, "fetch_semantic_scholar_batch", fake_fetch_semantic_scholar_batch
    )
    monkeypatch.setattr(pl, "fetch_openalex", lambda *a, **k: {"OpenAlex.Error": ""})
    monkeypatch.setattr(pl, "fetch_crossref", lambda *a, **k: {"crossref.Error": ""})

    def fake_get_limiter(name, rps, burst=None):
        calls["limiters"].append((name, rps, burst))
        return DummyLimiter()

    monkeypatch.setattr(pl, "get_limiter", fake_get_limiter)

    base_cfg = pl.Config()
    cfg_dict = base_cfg.model_dump()
    cfg_dict["system"]["rate"]["global_rps"] = 8
    cfg_dict["system"]["rate"]["global_burst"] = 9
    cfg_dict["sources"]["pubmed"]["rps"] = 5
    cfg_dict["sources"]["pubmed"]["burst"] = 6
    custom_cfg = pl.Config.model_validate(cfg_dict)

    def fake_apply_config_overrides(*args, **kwargs):
        return custom_cfg

    monkeypatch.setattr(pl.cli, "apply_config_overrides", fake_apply_config_overrides)

    input_csv = Path("tests/data/pmids.csv")
    output_csv = tmp_path / "override.csv"

    exit_code = pl.main(
        [
            "--input-csv",
            str(input_csv),
            "--output",
            str(output_csv),
            "--log-level",
            "ERROR",
        ]
    )

    assert exit_code == 0
    assert output_csv.exists()
    assert calls["limiters"][0] == ("pubmed", 5, 6)
    assert delays["pubmed"] == pytest.approx(1 / 5)
    assert delays["semantic_scholar"] == pytest.approx(1 / 5)
