import logging

import pytest

from library.clients import pubchem
from library.config import PubChemCfg
from library.integration import pubchem_library


@pytest.mark.unit
def test_resolve_pubchem_record__service_unavailable_stops_resolution(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    cfg = PubChemCfg()
    cfg.resolve_order = ("smiles", "inchikey")

    identifiers = {"canonical_smiles": "CCO", "standard_inchi_key": "ABC"}

    resolution_cache: dict[object, pubchem_library.PubChemResolution] = {}
    resolution_key = ("key",)

    call_log: list[str] = []

    def fail_smiles(value: str, cfg_arg: PubChemCfg) -> str | None:
        call_log.append("smiles")
        raise pubchem.PubChemServiceUnavailable(
            "server_error", {"status": 503, "retry_after": 30.0}
        )

    def should_not_call(*_args, **_kwargs) -> str | None:
        pytest.fail("Subsequent resolvers must not be invoked")

    monkeypatch.setattr(pubchem_library, "get_cid_from_smiles", fail_smiles)
    monkeypatch.setattr(pubchem_library, "get_cid_from_inchikey", should_not_call)

    caplog.set_level("WARNING", logger="chembl")

    resolution = pubchem_library.resolve_pubchem_record(
        identifiers,
        cfg,
        resolution_cache=resolution_cache,
        resolution_key=resolution_key,
    )

    assert resolution.cid is None
    assert resolution.temporary_failure is True
    assert resolution.status == 503
    assert call_log == ["smiles"]
    assert resolution_cache[resolution_key] == resolution
    assert any(
        record.message.startswith("pubchem_unavailable") for record in caplog.records
    )


@pytest.mark.unit
def test_resolve_pubchem_record__service_unavailable_cooldown_logs_info(
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    cfg = PubChemCfg()
    cfg.resolve_order = ("smiles",)

    identifiers = {"canonical_smiles": "CCO"}

    details = {
        "status": 503,
        "cooldown_remaining": 12.5,
        "cooldown_started_at": 100.0,
        "cooldown_until": 112.5,
        "retry_after": 30.0,
        "retry_after_source": "header",
        "cache": True,
    }

    def fail_with_cooldown(value: str, cfg_arg: PubChemCfg) -> str | None:
        raise pubchem.PubChemServiceUnavailable("server_error", details)

    monkeypatch.setattr(pubchem_library, "get_cid_from_smiles", fail_with_cooldown)

    caplog.set_level("INFO", logger="chembl")

    resolution = pubchem_library.resolve_pubchem_record(identifiers, cfg)

    assert resolution.cid is None
    assert resolution.temporary_failure is True
    assert resolution.status == 503
    assert any(
        record.levelno == logging.INFO
        and record.message.startswith("pubchem_unavailable")
        for record in caplog.records
    )
    assert not any(
        record.levelno >= logging.WARNING
        and record.message.startswith("pubchem_unavailable")
        for record in caplog.records
    )
