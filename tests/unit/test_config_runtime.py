"""Tests for :mod:`library.config.runtime` helpers."""

from __future__ import annotations

from library.config.models import ApiCfg, CrossRefCfg, OpenAlexCfg, RetryCfg
from library.config.runtime import crossref_session, openalex_session


def test_openalex_session__respects_verify_flag() -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg()
    cfg = OpenAlexCfg(mailto="chembl-testing@ebi.ac.uk", verify=False)

    with openalex_session(api_cfg, retry_cfg, cfg) as session:
        assert session.verify is False
        assert session.headers["mailto"] == cfg.mailto


def test_crossref_session__accepts_bundle_path(tmp_path) -> None:
    api_cfg = ApiCfg()
    retry_cfg = RetryCfg()
    bundle = tmp_path / "bundle.pem"
    bundle.write_text("dummy", encoding="utf-8")
    cfg = CrossRefCfg(mailto="chembl-testing@ebi.ac.uk", verify=bundle.as_posix())

    with crossref_session(api_cfg, retry_cfg, cfg) as session:
        assert session.verify == bundle.as_posix()
        assert session.headers["mailto"] == cfg.mailto
