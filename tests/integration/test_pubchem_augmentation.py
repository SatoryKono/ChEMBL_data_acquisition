from __future__ import annotations

import json
import logging
from collections.abc import Hashable, Mapping, MutableMapping
from pathlib import Path

import pandas as pd
import pytest

from library.config import PubChemCfg
from library.pipelines.testitem import cli as testitem_cli
from library.pipelines.testitem import pubchem
from library.pipelines.testitem.catalog import ParentLookupStats


@pytest.mark.integration
def test_add_pubchem_data__normal_flow_populates_cache(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"
    pubchem_cfg.cache_ttl_hours = None

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "molecule_type": ["Small molecule", "Small molecule"],
        }
    )

    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[str, object] = {}
    parent_record_cache: dict[str, pd.Series | None] = {}

    def fake_resolve(
        frame_arg: pd.DataFrame,
        cfg_arg: PubChemCfg,
        *,
        cid_cache: dict[str, str | None],
        resolution_cache: dict[str, object],
        load_parent_record,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> tuple[pd.Series, set[str], bool]:
        assert frame_arg.equals(frame)
        assert skip_mask.tolist() == [False, False]
        assert prefer_local_mask.tolist() == [False, False]
        cid_series = pd.Series(["CID1", "CID2"], index=frame.index, dtype="string")
        cid_cache["CHEMBL1"] = "CID1"
        cid_cache["CHEMBL2"] = "CID2"
        return cid_series, {"CID1", "CID2"}, True

    def fake_merge(
        frame_arg: pd.DataFrame,
        cid_series: pd.Series,
        lookup_cids: set[str],
        *,
        cfg: PubChemCfg,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> pd.DataFrame:
        assert lookup_cids == {"CID1", "CID2"}
        assert skip_mask.tolist() == [False, False]
        return pd.DataFrame(
            {
                "pubchem_cid": cid_series.astype("string"),
                "pubchem_iupac_name": ["Name1", "Name2"],
                "pubchem_canonical_smiles": ["SMILES1", "SMILES2"],
            },
            index=frame.index,
        )

    monkeypatch.setattr(pubchem, "_resolve_pubchem_cids", fake_resolve)
    monkeypatch.setattr(pubchem, "_merge_pubchem_properties", fake_merge)

    result = pubchem.add_pubchem_data(
        frame,
        pubchem_cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        parent_record_cache=parent_record_cache,
    )

    assert result["pubchem_cid"].tolist() == ["CID1", "CID2"]
    assert result["pubchem_iupac_name"].tolist() == ["Name1", "Name2"]
    assert cid_cache == {"CHEMBL1": "CID1", "CHEMBL2": "CID2"}
    assert pubchem_cfg.cid_cache_path is not None
    cache_payload = json.loads(pubchem_cfg.cid_cache_path.read_text(encoding="utf-8"))
    assert cache_payload["values"] == {"CHEMBL1": "CID1", "CHEMBL2": "CID2"}


@pytest.mark.integration
def test_add_pubchem_data__makes_pubchem_requests_when_enabled(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "canonical_smiles": ["CC(=O)OC1=CC=CC=C1C(=O)O"],
        }
    )

    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[str, object] = {}
    parent_record_cache: dict[str, pd.Series | None] = {}

    from library.integration import pubchem_library

    resolve_calls: list[dict[str, object]] = []

    def fake_resolve(
        identifiers: Mapping[str, str | None],
        cfg_arg: PubChemCfg,
        *,
        cid_cache: MutableMapping[str, str | None] | None = None,
        cache_key: str | None = None,
        resolution_cache: (
            MutableMapping[Hashable, pubchem_library.PubChemResolution] | None
        ) = None,
        resolution_key: Hashable | None = None,
    ) -> pubchem_library.PubChemResolution:
        resolve_calls.append({"identifiers": dict(identifiers), "cache_key": cache_key})
        assert cache_key == "CHEMBL1"
        assert identifiers.get("canonical_smiles") == "CC(=O)OC1=CC=CC=C1C(=O)O"
        if cid_cache is not None and cache_key is not None:
            cid_cache[cache_key] = "CID321"
        return pubchem_library.PubChemResolution(cid="CID321", source="smiles")

    properties_calls: list[str] = []

    def fake_get_properties(
        cid: str, cfg_arg: PubChemCfg
    ) -> pubchem_library.Properties:
        properties_calls.append(cid)
        assert cid == "CID321"
        return pubchem_library.Properties(
            "Acetylsalicylic acid",
            "C9H8O4",
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "InChI=1S/C9H8O4/...",
            "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        )

    monkeypatch.setattr(pubchem_library, "resolve_pubchem_record", fake_resolve)
    monkeypatch.setattr(pubchem_library, "get_properties", fake_get_properties)

    result = pubchem.add_pubchem_data(
        frame,
        pubchem_cfg,
        cid_cache=cid_cache,
        resolution_cache=resolution_cache,
        parent_record_cache=parent_record_cache,
    )

    assert resolve_calls, "Expected PubChem resolution to be triggered"
    assert properties_calls == ["CID321"]
    assert result.loc[0, "pubchem_cid"] == "CID321"
    assert result.loc[0, "pubchem_iupac_name"] == "Acetylsalicylic acid"
    assert result.loc[0, "pubchem_molecular_formula"] == "C9H8O4"
    assert cid_cache == {"CHEMBL1": "CID321"}


@pytest.mark.integration
def test_add_pubchem_data__skips_polymers_when_disallowed(
    tmp_path: Path,
    cfg,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.allow_polymer = False
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
            "molecule_type": ["Polymer", "Mixture", "Small molecule"],
        }
    )

    def fake_resolve(
        frame_arg: pd.DataFrame,
        cfg_arg: PubChemCfg,
        *,
        cid_cache: dict[str, str | None],
        resolution_cache: dict[str, object],
        load_parent_record,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> tuple[pd.Series, set[str], bool]:
        assert skip_mask.tolist() == [True, True, False]
        assert prefer_local_mask.tolist() == [False, False, False]
        cid_series = pd.Series(
            [pd.NA, pd.NA, "CID3"], index=frame.index, dtype="string"
        )
        cid_cache.update({})
        return cid_series, {"CID3"}, False

    def fake_merge(
        frame_arg: pd.DataFrame,
        cid_series: pd.Series,
        lookup_cids: set[str],
        *,
        cfg: PubChemCfg,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> pd.DataFrame:
        assert lookup_cids == {"CID3"}
        assert skip_mask.tolist() == [True, True, False]
        return pd.DataFrame(
            {
                "pubchem_cid": cid_series,
                "pubchem_iupac_name": [pd.NA, pd.NA, "Name3"],
            },
            index=frame.index,
        )

    monkeypatch.setattr(pubchem, "_resolve_pubchem_cids", fake_resolve)
    monkeypatch.setattr(pubchem, "_merge_pubchem_properties", fake_merge)

    events: list[tuple[str, dict[str, object]]] = []

    def capture_warning(event: str, **fields: object) -> None:
        logging.getLogger("pubchem-test").warning(event)
        events.append((event, fields))

    monkeypatch.setattr(pubchem.logger, "warning", capture_warning)
    caplog.set_level("WARNING", logger="pubchem-test")
    result = pubchem.add_pubchem_data(frame, pubchem_cfg)

    assert pd.isna(result.loc[0, "pubchem_cid"])
    assert pd.isna(result.loc[1, "pubchem_cid"])
    assert result.loc[2, "pubchem_cid"] == "CID3"
    assert any(event == "pubchem_skip_polymers" for event, _ in events)
    assert any(
        record.message == "pubchem_skip_polymers"
        for record in caplog.records
        if record.name == "pubchem-test"
    )


@pytest.mark.integration
def test_add_pubchem_data__uses_dictionary_mapping(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"

    dictionary_path = (
        Path(__file__).resolve().parents[1]
        / "resources"
        / "dictionaries"
        / "_testitem"
        / "chembl_to_pubchem.csv"
    )
    mapping_frame = pd.read_csv(dictionary_path)

    cid_cache: dict[str, str | None] = {
        str(row.molecule_chembl_id).strip().upper(): str(row.pubchem_cid).strip()
        for row in mapping_frame.itertuples()
    }
    expected_cache = cid_cache.copy()

    from library.integration import pubchem_library

    def forbid_resolve(*_args, **_kwargs) -> None:
        raise AssertionError("resolve_pubchem_cid should not be invoked for cached IDs")

    property_calls: list[str] = []

    def fake_get_properties(
        cid: str, cfg_arg: PubChemCfg
    ) -> pubchem_library.Properties:
        property_calls.append(cid)
        return pubchem_library.Properties(
            f"Name-{cid}",
            f"Formula-{cid}",
            f"iSMILES-{cid}",
            f"cSMILES-{cid}",
            f"InChI-{cid}",
            f"InChIKey-{cid}",
        )

    monkeypatch.setattr(pubchem, "resolve_pubchem_cid", forbid_resolve)
    monkeypatch.setattr(pubchem_library, "get_properties", fake_get_properties)

    frame = pd.DataFrame({"molecule_chembl_id": mapping_frame["molecule_chembl_id"]})

    result = pubchem.add_pubchem_data(
        frame,
        pubchem_cfg,
        cid_cache=cid_cache,
        resolution_cache={},
        parent_record_cache={},
    )

    expected_cids = mapping_frame["pubchem_cid"].tolist()
    assert result["pubchem_cid"].tolist() == expected_cids
    assert result["pubchem_iupac_name"].tolist() == [
        f"Name-{cid}" for cid in expected_cids
    ]
    assert result["pubchem_molecular_formula"].tolist() == [
        f"Formula-{cid}" for cid in expected_cids
    ]
    assert cid_cache == expected_cache
    assert set(property_calls) == set(expected_cids)


@pytest.mark.integration
def test_merge_pubchem_properties__preserves_existing_values_on_skip(cfg) -> None:
    pubchem_cfg = cfg.pubchem
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "pubchem_cid": ["CID1"],
            "pubchem_iupac_name": ["Existing"],
            "pubchem_canonical_smiles": ["SMILES"],
            "pubchem_isomeric_smiles": ["ISMILES"],
            "pubchem_inchi": ["InChI"],
            "pubchem_inchikey": ["InChIKey"],
            "pubchem_molecular_formula": ["C2H4O2"],
        }
    )

    skip_mask = pd.Series([True], index=frame.index, dtype="bool")
    prefer_local_mask = pd.Series([False], index=frame.index, dtype="bool")

    merged = pubchem._merge_pubchem_properties(
        frame,
        frame["pubchem_cid"].astype("string"),
        set(),
        cfg=pubchem_cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    for column in pubchem.PUBCHEM_COLUMNS:
        assert merged.loc[0, column] == frame.loc[0, column]


@pytest.mark.integration
def test_add_pubchem_data__keeps_existing_when_lookup_empty(
    cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "pubchem_cid": ["CID1"],
            "pubchem_canonical_smiles": ["LOCAL_SMILES"],
        }
    )

    def fake_prepare_caches(*_args, **_kwargs):
        return {}, {}, {}, set(), lambda *_a, **_kw: None

    def fake_prefetch(*_args, **_kwargs):
        return None

    def fake_resolve(*_args, **_kwargs):
        cid_series = pd.Series(["CID1"], index=frame.index, dtype="string")
        return cid_series, set(), False

    def fake_merge(
        frame_arg: pd.DataFrame,
        cid_series: pd.Series,
        lookup_cids: set[str],
        *,
        cfg,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> pd.DataFrame:
        assert lookup_cids == set()
        assert skip_mask.tolist() == [False]
        assert prefer_local_mask.tolist() == [False]
        return pd.DataFrame(
            {
                "pubchem_cid": cid_series,
                "pubchem_canonical_smiles": pd.Series(
                    [pd.NA], index=frame_arg.index, dtype="string"
                ),
            }
        )

    monkeypatch.setattr(pubchem, "_prepare_pubchem_caches", fake_prepare_caches)
    monkeypatch.setattr(pubchem, "_prefetch_parents", fake_prefetch)
    monkeypatch.setattr(pubchem, "_resolve_pubchem_cids", fake_resolve)
    monkeypatch.setattr(pubchem, "_merge_pubchem_properties", fake_merge)
    monkeypatch.setattr(
        pubchem, "_write_pubchem_cid_cache", lambda *args, **kwargs: None
    )

    result = pubchem.add_pubchem_data(frame, pubchem_cfg)

    assert result.loc[0, "pubchem_cid"] == "CID1"
    assert result.loc[0, "pubchem_canonical_smiles"] == "LOCAL_SMILES"


def test_merge_pubchem_properties__retains_partial_values_on_failed_lookup(
    cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.prefer_local_values = True

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "pubchem_cid": ["CID_OLD"],
            "pubchem_iupac_name": ["Existing"],
            "pubchem_canonical_smiles": ["SMILES_EXISTING"],
            "pubchem_isomeric_smiles": ["ISMILES_EXISTING"],
        }
    )

    cid_series = pd.Series(["CID_NEW"], index=frame.index, dtype="string")
    lookup_cids = {"CID_NEW"}

    from library.integration import pubchem_library

    class _DummyPubChemLib:
        Properties = pubchem_library.Properties

        def get_properties(
            self, cid: str, cfg_arg: PubChemCfg
        ) -> pubchem_library.Properties:
            assert cid == "CID_NEW"
            return pubchem_library.Properties(None, None, None, None, None, None)

    monkeypatch.setattr(pubchem, "_load_pubchem_library", lambda: _DummyPubChemLib())

    skip_mask = pd.Series([False], index=frame.index, dtype="bool")
    prefer_local_mask = pd.Series([False], index=frame.index, dtype="bool")

    merged = pubchem._merge_pubchem_properties(
        frame,
        cid_series,
        lookup_cids,
        cfg=pubchem_cfg,
        skip_mask=skip_mask,
        prefer_local_mask=prefer_local_mask,
    )

    assert merged.loc[0, "pubchem_cid"] == "CID_NEW"
    assert merged.loc[0, "pubchem_iupac_name"] == "Existing"
    assert merged.loc[0, "pubchem_canonical_smiles"] == "SMILES_EXISTING"
    assert merged.loc[0, "pubchem_isomeric_smiles"] == "ISMILES_EXISTING"


def test_merge_pubchem_properties__stops_after_service_unavailable(
    cfg, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.batch_size = 1
    pubchem_cfg.rps = 1

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
        }
    )

    cid_series = pd.Series(["CID1", "CID2", "CID3"], index=frame.index, dtype="string")
    lookup_cids = set(cid_series.tolist())

    skip_mask = pd.Series([False, False, False], index=frame.index, dtype="bool")
    prefer_local_mask = pd.Series(
        [False, False, False], index=frame.index, dtype="bool"
    )

    from library.integration import pubchem_library

    call_count = 0

    class _DummyPubChemLib:
        Properties = pubchem_library.Properties
        PubChemServiceUnavailable = pubchem_library.PubChemServiceUnavailable

        def get_properties(
            self, cid: str, cfg_arg: PubChemCfg
        ) -> pubchem_library.Properties:
            nonlocal call_count
            call_count += 1
            raise pubchem_library.PubChemServiceUnavailable(
                "server_error", {"status": 503, "retry_after": 30.0}
            )

    monkeypatch.setattr(pubchem, "_load_pubchem_library", lambda: _DummyPubChemLib())

    with caplog.at_level(logging.WARNING):
        merged = pubchem._merge_pubchem_properties(
            frame,
            cid_series,
            lookup_cids,
            cfg=pubchem_cfg,
            skip_mask=skip_mask,
            prefer_local_mask=prefer_local_mask,
        )

    assert call_count == 1
    assert merged["pubchem_cid"].tolist() == ["CID1", "CID2", "CID3"]
    remaining_columns = [col for col in pubchem.PUBCHEM_COLUMNS if col != "pubchem_cid"]
    assert merged[remaining_columns].isna().all().all()

    warning_messages = [
        record.message for record in caplog.records if record.levelno >= logging.WARNING
    ]
    failed_messages = [
        msg for msg in warning_messages if msg.startswith("pubchem_properties_failed")
    ]
    assert len(failed_messages) == 1
    assert any("retry_after=30.0" in msg for msg in failed_messages)
    assert any("status=503" in msg for msg in failed_messages)
    unavailable_messages = [
        msg
        for msg in warning_messages
        if msg.startswith("pubchem_properties_unavailable")
    ]
    assert unavailable_messages
    assert any("pending=3" in msg for msg in unavailable_messages)
    assert any("retry_after=30.0" in msg for msg in unavailable_messages)
    assert any("outcome='server_error'" in msg for msg in unavailable_messages)


def test_merge_pubchem_properties__limits_outstanding_requests_on_service_unavailable(
    cfg, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.batch_size = 10
    pubchem_cfg.rps = 3

    frame = pd.DataFrame(
        {
            "molecule_chembl_id": [f"CHEMBL{i}" for i in range(1, 8)],
        }
    )

    cid_series = pd.Series(
        [f"CID{i}" for i in range(1, 8)], index=frame.index, dtype="string"
    )
    lookup_cids = set(cid_series.tolist())

    skip_mask = pd.Series([False] * len(frame), index=frame.index, dtype="bool")
    prefer_local_mask = pd.Series([False] * len(frame), index=frame.index, dtype="bool")

    from library.integration import pubchem_library

    call_count = 0

    class _DummyPubChemLib:
        Properties = pubchem_library.Properties
        PubChemServiceUnavailable = pubchem_library.PubChemServiceUnavailable

        def get_properties(
            self, cid: str, cfg_arg: PubChemCfg
        ) -> pubchem_library.Properties:
            nonlocal call_count
            call_count += 1
            raise pubchem_library.PubChemServiceUnavailable(
                "server_error", {"status": 503, "retry_after": 30.0}
            )

    monkeypatch.setattr(pubchem, "_load_pubchem_library", lambda: _DummyPubChemLib())

    with caplog.at_level(logging.WARNING):
        merged = pubchem._merge_pubchem_properties(
            frame,
            cid_series,
            lookup_cids,
            cfg=pubchem_cfg,
            skip_mask=skip_mask,
            prefer_local_mask=prefer_local_mask,
        )

    assert call_count == pubchem_cfg.rps
    assert merged["pubchem_cid"].tolist() == [f"CID{i}" for i in range(1, 8)]
    remaining_columns = [col for col in pubchem.PUBCHEM_COLUMNS if col != "pubchem_cid"]
    assert merged[remaining_columns].isna().all().all()

    warning_messages = [
        record.message for record in caplog.records if record.levelno >= logging.WARNING
    ]
    assert (
        sum(msg.startswith("pubchem_properties_failed") for msg in warning_messages)
        == 1
    )
    unavailable_messages = [
        msg
        for msg in warning_messages
        if msg.startswith("pubchem_properties_unavailable")
    ]
    assert unavailable_messages
    assert any("pending=7" in msg for msg in unavailable_messages)


@pytest.mark.integration
def test_resolve_pubchem_cid__temporary_failure_does_not_cache(
    cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    from library.integration import pubchem_library

    row = pd.Series({"molecule_chembl_id": "CHEMBL1"})
    cid_cache: dict[str, str | None] = {}
    resolution_cache: dict[str, object] = {}

    class _DummyPubChemLib:
        def resolve_pubchem_record(self, *args, **kwargs):
            return pubchem_library.PubChemResolution(
                cid=None,
                source=None,
                status=503,
                temporary_failure=True,
            )

    loader_called = False

    def fake_parent_loader(_identifier: str) -> None:
        nonlocal loader_called
        loader_called = True
        return None

    monkeypatch.setattr(pubchem, "_load_pubchem_library", lambda: _DummyPubChemLib())

    result = pubchem.resolve_pubchem_cid(
        row,
        cid_cache,
        pubchem_cfg,
        parent_loader=fake_parent_loader,
        resolution_cache=resolution_cache,
    )

    assert result is None
    assert cid_cache == {}
    assert loader_called is False


@pytest.mark.integration
def test_augment_pubchem__initialises_session_and_reuses_cache(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.cache_ttl_hours = None
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"

    initial_payload = {
        "metadata": {"version": 1, "updated_at": "2020-01-01T00:00:00+00:00"},
        "values": {"CHEMBL1": "CID1"},
    }
    pubchem_cfg.cid_cache_path.write_text(
        json.dumps(initial_payload, indent=2), encoding="utf-8"
    )

    init_calls: list[tuple[object, object]] = []

    def record_init(api_cfg: object, retry_cfg: object) -> None:
        init_calls.append((api_cfg, retry_cfg))

    monkeypatch.setattr(pubchem.pl, "init_session", record_init)
    monkeypatch.setattr(pubchem, "_PUBCHEM_SESSION_SIGNATURE", None)

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1", "CHEMBL2"]}, index=[0, 1])

    call_count = 0

    def fake_resolve(
        frame_arg: pd.DataFrame,
        cfg_arg: PubChemCfg,
        *,
        cid_cache: dict[str, str | None],
        resolution_cache: dict[str, object],
        load_parent_record,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> tuple[pd.Series, set[str], bool]:
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            assert cid_cache.get("CHEMBL1") == "CID1"
            assert "CHEMBL2" not in cid_cache
            cid_cache["CHEMBL2"] = "CID2"
            series = pd.Series(["CID1", "CID2"], index=frame.index, dtype="string")
            return series, {"CID2"}, True
        assert cid_cache.get("CHEMBL1") == "CID1"
        assert cid_cache.get("CHEMBL2") == "CID2"
        series = pd.Series(["CID1", "CID2"], index=frame.index, dtype="string")
        return series, set(), False

    def fake_merge(
        frame_arg: pd.DataFrame,
        cid_series: pd.Series,
        lookup_cids: set[str],
        *,
        cfg: PubChemCfg,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "pubchem_cid": cid_series,
                "pubchem_iupac_name": ["Name1", "Name2"],
            },
            index=frame.index,
        )

    monkeypatch.setattr(pubchem, "_resolve_pubchem_cids", fake_resolve)
    monkeypatch.setattr(pubchem, "_merge_pubchem_properties", fake_merge)

    client = object()
    first = pubchem.augment_pubchem(
        frame,
        pubchem_cfg=pubchem_cfg,
        api_cfg=cfg.api,
        retry_cfg=cfg.system.retry,
        timeout=5.0,
        client=client,
        fields=None,
        request_limit=10,
    )
    second = pubchem.augment_pubchem(
        frame,
        pubchem_cfg=pubchem_cfg,
        api_cfg=cfg.api,
        retry_cfg=cfg.system.retry,
        timeout=5.0,
        client=client,
        fields=None,
        request_limit=10,
    )

    assert first["pubchem_cid"].tolist() == ["CID1", "CID2"]
    assert second["pubchem_cid"].tolist() == ["CID1", "CID2"]
    assert len(init_calls) == 1
    cache_payload = json.loads(pubchem_cfg.cid_cache_path.read_text(encoding="utf-8"))
    assert cache_payload["values"] == {"CHEMBL1": "CID1", "CHEMBL2": "CID2"}


@pytest.mark.integration
@pytest.mark.pipeline_scenario("enrichment", "missing_data")
def test_finalize_output__pubchem_fallback_partial_rows(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg.system.doc_quality.enable = False

    chunk = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(
                ["CHEMBL_OK", "CHEMBL_MISSING"], dtype="string"
            ),
            "natural_product": pd.Series([True, False], dtype="boolean"),
            "pubchem_cid": pd.Series(["CID_OK", pd.NA], dtype="object"),
            "pubchem_iupac_name": pd.Series(
                ["Existing", pd.NA], dtype="string"
            ),
        }
    )

    def _empty_report(*_: object, **__: object) -> pd.DataFrame:
        return pd.DataFrame({"column": pd.Series(dtype="string")})

    monkeypatch.setattr(
        testitem_cli.qc_report,
        "generate_qc_report",
        _empty_report,
    )
    monkeypatch.setattr(
        testitem_cli.data_correlation,
        "generate_correlation_report",
        _empty_report,
    )

    augment_calls: list[pd.DataFrame] = []

    def _fake_pubchem_augment(frame: pd.DataFrame, **_: object) -> pd.DataFrame:
        augment_calls.append(frame.copy())
        augmented = frame.copy()
        augmented.loc[:, "pubchem_cid"] = pd.Series(
            ["CID_FALLBACK"], index=frame.index, dtype="string"
        )
        augmented.loc[:, "pubchem_iupac_name"] = pd.Series(
            ["Recovered"], index=frame.index, dtype="string"
        )
        return augmented

    monkeypatch.setattr(
        testitem_cli,
        "_load_pubchem_augmenter",
        lambda: _fake_pubchem_augment,
    )

    class _StubChemblClient:
        pass

    pubchem_context = testitem_cli.PubChemAugmentationContext(
        pubchem_cfg=cfg.pubchem,
        api_cfg=cfg.api,
        retry_cfg=cfg.retry,
        client=_StubChemblClient(),
        timeout=5.0,
        fields=("pubchem_cid", "pubchem_iupac_name"),
        request_limit=10,
    )

    stats = ParentLookupStats(
        source="lookup",
        missing=0,
        unique=2,
        attached=2,
        uncovered=0,
        hierarchy_attached=0,
        fallback_attached=0,
        no_parent=0,
    )

    def _stats_supplier() -> ParentLookupStats:
        return stats

    output_path = tmp_path / "partial_fallback.csv"

    exit_code, artifacts = testitem_cli.finalize_output(
        [chunk],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=_stats_supplier,
        input_csv=output_path,
        pubchem_context=pubchem_context,
    )

    assert exit_code == 0
    assert artifacts is not None
    assert len(augment_calls) == 1

    augmented_frame = augment_calls[0]
    assert list(augmented_frame["molecule_chembl_id"]) == ["CHEMBL_MISSING"]
    assert augmented_frame.index.tolist() == [1]

    dataset = pd.read_csv(
        artifacts.dataset,
        dtype={
            "molecule_chembl_id": "string",
            "pubchem_cid": "string",
            "pubchem_iupac_name": "string",
        },
    )
    dataset = dataset.set_index("molecule_chembl_id", drop=True)

    assert dataset.loc["CHEMBL_OK", "pubchem_cid"] == "CID_OK"
    assert dataset.loc["CHEMBL_OK", "pubchem_iupac_name"] == "Existing"
    assert dataset.loc["CHEMBL_MISSING", "pubchem_cid"] == "CID_FALLBACK"
    assert dataset.loc["CHEMBL_MISSING", "pubchem_iupac_name"] == "Recovered"


@pytest.mark.integration
def test_add_pubchem_data__expires_stale_cache(
    tmp_path: Path,
    cfg,
    monkeypatch: pytest.MonkeyPatch,
    caplog: pytest.LogCaptureFixture,
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = True
    pubchem_cfg.cache_ttl_hours = 1
    pubchem_cfg.cid_cache_path = tmp_path / "cid_cache.json"

    stale_payload = {
        "metadata": {"version": 1, "updated_at": "2019-12-30T00:00:00+00:00"},
        "values": {"CHEMBL1": "CID1"},
    }
    pubchem_cfg.cid_cache_path.write_text(
        json.dumps(stale_payload, indent=2), encoding="utf-8"
    )

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL2"]})

    def fake_resolve(
        frame_arg: pd.DataFrame,
        cfg_arg: PubChemCfg,
        *,
        cid_cache: dict[str, str | None],
        resolution_cache: dict[str, object],
        load_parent_record,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> tuple[pd.Series, set[str], bool]:
        assert cid_cache == {}
        series = pd.Series(["CID2"], index=frame.index, dtype="string")
        cid_cache["CHEMBL2"] = "CID2"
        return series, {"CID2"}, True

    def fake_merge(
        frame_arg: pd.DataFrame,
        cid_series: pd.Series,
        lookup_cids: set[str],
        *,
        cfg: PubChemCfg,
        skip_mask: pd.Series,
        prefer_local_mask: pd.Series,
    ) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "pubchem_cid": cid_series,
                "pubchem_iupac_name": ["Name2"],
            },
            index=frame.index,
        )

    monkeypatch.setattr(pubchem, "_resolve_pubchem_cids", fake_resolve)
    monkeypatch.setattr(pubchem, "_merge_pubchem_properties", fake_merge)

    info_events: list[str] = []

    def capture_info(event: str, **fields: object) -> None:
        logging.getLogger("pubchem-test").info(event)
        info_events.append(event)

    monkeypatch.setattr(pubchem.logger, "info", capture_info)
    caplog.set_level("INFO", logger="pubchem-test")
    result = pubchem.add_pubchem_data(frame, pubchem_cfg)

    assert result["pubchem_cid"].tolist() == ["CID2"]
    assert "pubchem_cache_expired" in info_events
    assert any(
        record.message == "pubchem_cache_expired"
        for record in caplog.records
        if record.name == "pubchem-test"
    )


@pytest.mark.integration
def test_add_pubchem_data__disabled_mode_passthrough(
    cfg, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    pubchem_cfg = cfg.pubchem
    pubchem_cfg.enable = False

    frame = pd.DataFrame({"molecule_chembl_id": ["CHEMBL1"]})

    info_events: list[str] = []

    def capture_info(event: str, **fields: object) -> None:
        logging.getLogger("pubchem-test").info(event)
        info_events.append(event)

    monkeypatch.setattr(pubchem.logger, "info", capture_info)
    caplog.set_level("INFO", logger="pubchem-test")
    result = pubchem.add_pubchem_data(frame, pubchem_cfg)

    assert result is frame
    assert "pubchem_disabled" in info_events
    for column in pubchem.PUBCHEM_COLUMNS:
        assert column in result.columns
        assert result[column].isna().all()
    assert any(
        record.message == "pubchem_disabled"
        for record in caplog.records
        if record.name == "pubchem-test"
    )
