"""Unit tests for the isoform post-processing helpers."""

from __future__ import annotations

from pathlib import Path
import warnings

import pandas as pd
import pytest

from library.postprocessing.target import isoform
from scripts import get_target_data


@pytest.mark.unit
def test_transform__supports_alias_columns():
    """Isoform transform accepts common legacy aliases for required columns."""

    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["SynA|SynB"],
            "isoform_names": ["Isoform 1|Isoform 2"],
            "isoform_ids": ["P12345-1|P12345-2"],
            "primary_accession": ["P12345"],
            "chembl_id": ["CHEMBL123"],
        }
    )

    result = isoform._transform(frame).result

    expected = pd.DataFrame(
        [
            {
                "id": "P12345-1",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "isoform 1",
            },
            {
                "id": "P12345-1",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "syna",
            },
            {
                "id": "P12345-2",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "isoform 2",
            },
            {
                "id": "P12345-2",
                "uniprot_id_primary": "P12345",
                "target_chembl_id": "CHEMBL123",
                "name": "synb",
            },
        ]
    )

    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.unit
def test_transform__missing_identifier_columns_warns_and_returns_empty_result():
    """A warning is emitted and an empty result returned when IDs are absent."""

    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["SynA"],
            "isoform_names": ["Isoform 1"],
            "isoform_ids": ["P12345-1"],
        }
    )

    with pytest.warns(UserWarning) as captured:
        result = isoform._transform(frame)

    assert len(captured) == 1
    message = str(captured[0].message)
    assert "uniprot_id_primary" in message
    assert "target_chembl_id" in message

    for attribute in (
        "result",
        "combined",
        "dedup_stage1",
        "sorted_stage",
        "dedup_stage2",
    ):
        frame_result = getattr(result, attribute)
        assert isinstance(frame_result, pd.DataFrame)
        assert frame_result.empty


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        "out.csv",
        "out_chembl.csv",
        "out_normalized.csv",
        "out_uniprot.csv",
        "output.targets.csv",
        "output.targets_normalized.csv",
    ],
)
def test_matches_expected_input_name__supports_modern_exports(filename: str) -> None:
    """Modern aggregated exports are accepted for isoform processing."""

    assert isoform._matches_expected_input_name(filename)


@pytest.mark.unit
@pytest.mark.parametrize(
    "filename",
    [
        ".output.targets_20251005.csv.tmp",
        "output.targets_20251005.csv.tmp",
        "output.targets_20251005.csv_normalized",
    ],
)
def test_matches_expected_input_name__normalises_temporary_variants(
    filename: str,
) -> None:
    """Temporary artefacts created during exports are treated as canonical."""

    assert isoform._matches_expected_input_name(filename)


@pytest.mark.unit
def test_resolve_input_path__temporary_export_emits_no_warning(tmp_path: Path) -> None:
    """Temporary exports with compound suffixes are accepted without warnings."""

    input_path = tmp_path / ".output.targets_20251005.csv.tmp"
    input_path.write_text("target_chembl_id\n", encoding="utf-8")

    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        resolved = isoform._resolve_input_path(input_path)

    assert resolved == input_path
    assert not captured


@pytest.mark.unit
def test_is_supported_target_export__accepts_cli_default(tmp_path: Path) -> None:
    """The CLI default ``output.targets_*.csv`` export is supported."""

    path = tmp_path / "output.targets_20250101.csv"
    path.write_text("target_chembl_id\n", encoding="utf-8")

    assert get_target_data._is_supported_target_export(path)


@pytest.mark.unit
def test_postprocess_target_exports__runs_helpers_for_cli_default(
    tmp_path: Path, cfg, monkeypatch: pytest.MonkeyPatch
) -> None:
    """All target post-processing helpers run for the CLI default export."""

    source = tmp_path / "output.targets_20250101.csv"
    source.write_text("target_chembl_id\n", encoding="utf-8")

    calls: list[str] = []

    def _organism_stub(source: Path, *, cfg) -> None:  # pragma: no cover - simple stub
        calls.append("organism")
        return None

    def _isoform_stub(
        source: Path,
        *,
        cfg,
        context=None,
        ambiguous_classifications=None,
    ) -> None:
        calls.append("isoform")
        return None

    def _names_stub(source: Path, *, cfg) -> None:  # pragma: no cover - simple stub
        calls.append("names")
        return None

    def _iuphar_stub(
        source: Path, *, verbose: bool = True
    ) -> None:  # pragma: no cover - stub
        calls.append("iuphar")
        return None

    monkeypatch.setattr(get_target_data, "_postprocess_organism_export", _organism_stub)
    monkeypatch.setattr(get_target_data, "_postprocess_isoform_export", _isoform_stub)
    monkeypatch.setattr(get_target_data, "_postprocess_names_export", _names_stub)
    monkeypatch.setattr(get_target_data, "_postprocess_iuphar_export", _iuphar_stub)

    get_target_data._postprocess_target_exports(
        source,
        cfg=cfg,
        context=None,
        ambiguous_classifications=None,
        verbose=True,
    )

    assert calls == ["organism", "isoform", "names", "iuphar"]
