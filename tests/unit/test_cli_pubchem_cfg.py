import pytest

from library.config import Config
from library.pipelines.testitem import cli


@pytest.mark.unit
def test_prepare_pubchem_api_cfg__uses_pubchem_contact(cfg: Config) -> None:
    api_cfg = cfg.api
    cfg.pubchem.user_agent = "owner@example.org"
    updated = cli._prepare_pubchem_api_cfg(cfg, api_cfg)
    assert updated.user_agent == "owner@example.org"


@pytest.mark.unit
def test_prepare_pubchem_api_cfg__raises_for_placeholder() -> None:
    cfg = Config()
    cfg.api.user_agent = "contact@example.org"
    cfg.pubchem.user_agent = "contact@example.org"
    with pytest.raises(ValueError):
        cli._prepare_pubchem_api_cfg(cfg, cfg.api)
