import argparse

import pytest

from library.cli import apply_config_overrides
from library.config import DEFAULT_CONFIG_PATH, Config, load_config


@pytest.mark.unit
def test_pubchem_enable__default_model_true() -> None:
    cfg = Config()
    assert cfg.pubchem.enable is True


@pytest.mark.unit
def test_pubchem_enable__config_file_true() -> None:
    cfg = load_config(DEFAULT_CONFIG_PATH)
    assert cfg.pubchem.enable is True


@pytest.mark.unit
def test_pubchem_enable__cli_defaults_true() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG_PATH)
    parser.add_argument("--timeout", type=float, default=30.0)
    parser.add_argument("--column", default="molecule_chembl_id")
    parser.add_argument("--batch_size", type=int, default=1000)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--offset", type=int, default=0)
    args = parser.parse_args([])

    mapping = {
        "timeout": "testitem.timeout",
        "column": "testitem.column",
        "batch_size": "testitem.batch_size",
        "limit": "testitem.limit",
        "offset": "testitem.offset",
    }

    cfg = apply_config_overrides(args, parser, args.config, mapping=mapping)

    assert cfg.pubchem.enable is True


@pytest.mark.unit
def test_pubchem_enable__env_alias_short(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("CHEMBL_DA_PUBCHEM_ENABLE", "false")

    cfg = load_config(DEFAULT_CONFIG_PATH)

    assert cfg.pubchem.enable is False
