from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

try:  # pragma: no cover - compatibility for pandera namespace changes
    import pandera.pandas as _pa  # noqa: F401
except ModuleNotFoundError:  # pragma: no cover - fallback for legacy pandera
    import pandera as _pa  # type: ignore[no-redef]

    sys.modules.setdefault("pandera.pandas", _pa)


import pandas as pd
import pytest

from library.pipelines.testitem import cli as testitem_cli
from library.pipelines.testitem.catalog import (
    PARENT_LOOKUP_SOURCE_SKIPPED,
    ParentLookupStats,
)

if TYPE_CHECKING:
    from library.config import Config


def _stub_parent_stats() -> ParentLookupStats:
    """Return a neutral ``ParentLookupStats`` snapshot for tests."""

    return ParentLookupStats(
        source=PARENT_LOOKUP_SOURCE_SKIPPED,
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    "output, fallback_date, expected",
    [
        (Path(".output.testitem_20240101.csv.tmp"), "19991231", ("testitem", "20240101")),
        ("output.testitems_20240101.csv", "19991231", ("testitem", "20240101")),
        ("output.testitem.csv", "19991231", ("testitem", "19991231")),
    ],
    )
def test_normalise_output_labels__derives_table_and_date(
    output: Path | str, fallback_date: str, expected: tuple[str, str]
) -> None:
    """Intermediate filenames should map to canonical artefact labels."""

    assert (
        testitem_cli._normalise_output_labels(output, fallback_date=fallback_date)  # noqa: SLF001
        == expected
    )


@pytest.mark.unit
@pytest.mark.pipeline_scenario("assembly")
def test_finalize_output__preserves_existing_optional_columns(
    tmp_path: Path, cfg: "Config", monkeypatch: pytest.MonkeyPatch
) -> None:
    """Existing optional column values must survive column alignment."""

    chunk_with_pubchem = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL1"],
            "pubchem_cid": ["123"],
        }
    )
    chunk_without_pubchem = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL2"],
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

    output_path = tmp_path / "output.csv"

    exit_code, artifacts = testitem_cli.finalize_output(
        [chunk_with_pubchem, chunk_without_pubchem],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=_stub_parent_stats,
        input_csv=output_path,
    )

    assert exit_code == 0
    assert artifacts is not None
    assert artifacts.dataset.exists()

    result_df = pd.read_csv(
        artifacts.dataset,
        dtype={"molecule_chembl_id": "string", "pubchem_cid": "string"},
    )
    cid_by_id = (
        result_df.set_index("molecule_chembl_id", drop=True)["pubchem_cid"].to_dict()
    )

    assert cid_by_id["CHEMBL1"] == "123"
    assert pd.isna(cid_by_id["CHEMBL2"])


@pytest.mark.unit
@pytest.mark.pipeline_scenario("export")
def test_finalize_output__applies_pubchem_fallback_after_postprocessing(
    tmp_path: Path, cfg: "Config", monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fallback PubChem enrichment should persist in the final export."""

    chunk_without_pubchem = pd.DataFrame(
        {
            "molecule_chembl_id": ["CHEMBL3"],
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
        augmented["pubchem_cid"] = pd.Series(["CID123"], dtype="string")
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
        timeout=10.0,
        fields=("pubchem_cid",),
        request_limit=50,
    )

    output_path = tmp_path / "fallback.csv"

    exit_code, artifacts = testitem_cli.finalize_output(
        [chunk_without_pubchem],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=_stub_parent_stats,
        input_csv=output_path,
        pubchem_context=pubchem_context,
    )

    assert exit_code == 0
    assert artifacts is not None
    assert len(augment_calls) == 1

    result_df = pd.read_csv(
        artifacts.dataset,
        dtype={"molecule_chembl_id": "string", "pubchem_cid": "string"},
    )

    cid_series = result_df.set_index("molecule_chembl_id", drop=True)["pubchem_cid"]
    assert cid_series.loc["CHEMBL3"] == "CID123"


@pytest.mark.unit
@pytest.mark.pipeline_scenario("export")
def test_finalize_output__triggers_pubchem_fallback_when_empty_strings(
    tmp_path: Path, cfg: "Config", monkeypatch: pytest.MonkeyPatch
) -> None:
    """Fallback PubChem enrichment must treat empty strings as missing values."""

    chunk_with_empty_pubchem = pd.DataFrame(
        {
            "molecule_chembl_id": pd.Series(["CHEMBL4"], dtype="string"),
            "pubchem_cid": pd.Series([""], dtype="string"),
            "pubchem_iupac_name": pd.Series([""], dtype="string"),
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
        augmented["pubchem_cid"] = pd.Series(["CID456"], dtype="string")
        augmented["pubchem_iupac_name"] = pd.Series(
            ["Acetylsalicylic acid"], dtype="string"
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
        timeout=10.0,
        fields=("pubchem_cid", "pubchem_iupac_name"),
        request_limit=50,
    )

    output_path = tmp_path / "fallback_empty.csv"

    exit_code, artifacts = testitem_cli.finalize_output(
        [chunk_with_empty_pubchem],
        cfg=cfg,
        output=output_path,
        parent_stats_supplier=_stub_parent_stats,
        input_csv=output_path,
        pubchem_context=pubchem_context,
    )

    assert exit_code == 0
    assert artifacts is not None
    assert len(augment_calls) == 1

    before_augment = augment_calls[0]
    assert pd.isna(before_augment.loc[0, "pubchem_cid"])
    assert pd.isna(before_augment.loc[0, "pubchem_iupac_name"])

    result_df = pd.read_csv(
        artifacts.dataset,
        dtype={
            "molecule_chembl_id": "string",
            "pubchem_cid": "string",
            "pubchem_iupac_name": "string",
        },
    )

    cid_series = result_df.set_index("molecule_chembl_id", drop=True)["pubchem_cid"]
    iupac_series = result_df.set_index("molecule_chembl_id", drop=True)[
        "pubchem_iupac_name"
    ]

    assert cid_series.loc["CHEMBL4"] == "CID456"
    assert iupac_series.loc["CHEMBL4"] == "Acetylsalicylic acid"
