import pytest

from library.config.models import OpenAlexCfg


@pytest.mark.unit
def test_openalex_cfg__defaults_limit_burst() -> None:
    cfg = OpenAlexCfg()

    assert cfg.rps == 4
    assert cfg.burst == 1
