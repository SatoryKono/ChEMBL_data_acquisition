"""Comprehensive tests for :mod:`library.postprocessing.activity_extended`."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import ActivityExtendedError, process_activity_extended
from library.postprocessing.activity_extended import (
    _FINAL_COLUMN_ORDER,
    _TARGET_COLUMNS,
    _apply_multimol_logic,
    _augment_activity_frame,
    _derive_output_path,
    _latest_activity_export,
    _load_citation_fraction,
    _load_document_lookup,
    _load_target_metadata,
    _prepare_unknown_chirality,
    _rename_columns,
    _resolve_targets_path,
    _select_and_cast,
    _transform_activity_frame,
)

pytestmark = pytest.mark.integration


EXPECTED_DTYPE_MAP: dict[str, str] = {
    "activity_chembl_id": "object",
    "saltform_id": "object",
    "molecule_chembl_id": "object",
    "molecule_chembl_id.1": "string",
    "target_chembl_id": "string",
    "assay_chembl_id": "object",
    "document_chembl_id": "object",
    "completed": "string",
    "bao_endpoint": "object",
    "standard_type": "object",
    "standard_value": "string",
    "pA_value": "string",
    "bao_format": "object",
    "compound_key": "object",
    "compound_name": "object",
    "standard_inchi_skeleton": "string",
    "unknown_chirality": "boolean",
    "multmol_assay": "boolean",
    "assay_with_same_target": "Int64",
    "exact_data_citation": "boolean",
    "higly_correlated_assay": "boolean",
    "shuffled_assay": "boolean",
    "review": "boolean",
    "rounded_data_citation": "boolean",
    "allosteric": "boolean",
    "nam": "boolean",
    "pam": "boolean",
    "high_citation_rate": "boolean",
    "original_activity_approx": "string",
    "original_activity_exact": "string",
    "is_citation": "boolean",
    "IUPHAR_class": "string",
    "IUPHAR_subclass": "string",
    "unicellular_organism": "boolean",
    "multifunctional_enzyme": "boolean",
    "gene_index": "string",
    "taxon_index": "string",
    "sortorder.target": "string",
}


@pytest.fixture()
def activity_resources(snapshot_resource: Path) -> Path:
    """Return the fixture directory with activity post-processing snapshots."""

    return snapshot_resource / "activity_extended"


def _copytree(src: Path, dst: Path) -> Path:
    shutil.copytree(src, dst)
    return dst


def test_load_document_lookup__renames_legacy_identifier(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    dictionary_root = tmp_path
    document_dir = dictionary_root / "_document"
    document_dir.mkdir(parents=True)

    frame = pd.DataFrame(
        {
            "ChEMBL.document_chembl_id": pd.Series(["DOC-1", "DOC-2"], dtype="string"),
            "completed": pd.Series(["yes", "no"], dtype="string"),
        }
    )
    frame.to_csv(
        document_dir / "document.csv", index=False, sep="\t", encoding="cp1252"
    )

    with caplog.at_level("WARNING"):
        lookup = _load_document_lookup(dictionary_root)

    assert lookup.columns.tolist() == ["document_chembl_id", "completed", "review"]
    assert lookup["document_chembl_id"].tolist() == ["DOC-1", "DOC-2"]
    assert lookup["review"].isna().all()
    assert "Renamed document identifier column" in caplog.text


def test_load_document_lookup__latin1_encoding_fallback(tmp_path: Path) -> None:
    dictionary_root = tmp_path
    document_dir = dictionary_root / "_document"
    document_dir.mkdir(parents=True)

    frame = pd.DataFrame(
        {
            "document_chembl_id": pd.Series(["DOC-1"], dtype="string"),
            "completed": pd.Series(["value\x81with_control"], dtype="string"),
        }
    )
    frame.to_csv(
        document_dir / "document.csv", index=False, sep="\t", encoding="latin-1"
    )

    lookup = _load_document_lookup(dictionary_root)

    assert lookup.loc[0, "document_chembl_id"] == "DOC-1"
    assert lookup.loc[0, "completed"] == "value\x81with_control"


def test_load_document_lookup__missing_identifier_reports_available_columns(
    tmp_path: Path,
) -> None:
    dictionary_root = tmp_path
    document_dir = dictionary_root / "_document"
    document_dir.mkdir(parents=True)

    pd.DataFrame({"identifier": ["DOC-1"]}).to_csv(
        document_dir / "document.csv", index=False, sep=",", encoding="utf-8"
    )

    with pytest.raises(ActivityExtendedError) as excinfo:
        _load_document_lookup(dictionary_root)

    message = str(excinfo.value)
    assert "identifier" in message
    assert "available columns" in message


def test_augment_activity_frame__fills_existing_blanks() -> None:
    frame = pd.DataFrame(
        {
            "activity_chembl_id": pd.Series(["", pd.NA], dtype="string"),
            "activity_id": pd.Series(["ACT-1", "ACT-2"], dtype="string"),
            "compound_name": pd.Series([pd.NA, ""], dtype="string"),
            "molecule_pref_name": pd.Series(
                ["Preferred-1", "Preferred-2"], dtype="string"
            ),
            "compound_key": pd.Series(["   ", pd.NA], dtype="string"),
            "molecule_chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
            "pchembl_value": pd.Series([5.0, 6.0], dtype="Float64"),
            "log_value": pd.Series([pd.NA, ""], dtype="object"),
        }
    )

    augmented, filled = _augment_activity_frame(frame)

    assert augmented["activity_chembl_id"].tolist() == ["ACT-1", "ACT-2"]
    assert augmented["compound_name"].tolist() == ["Preferred-1", "Preferred-2"]
    assert augmented["compound_key"].tolist() == ["CHEMBL1", "CHEMBL2"]
    pd.testing.assert_series_equal(
        augmented["log_value"],
        pd.Series([5.0, 6.0], index=frame.index, dtype="Float64"),
        check_names=False,
    )
    assert {
        "activity_chembl_id",
        "compound_name",
        "compound_key",
        "log_value",
    } <= filled


def test_augment_activity_frame__fills_blank_salt_and_log_value() -> None:
    frame = pd.DataFrame(
        {
            "salt_chembl_id": pd.Series(["", pd.NA], dtype="string"),
            "molecule_chembl_id": pd.Series(["CHEMBL1", "CHEMBL2"], dtype="string"),
            "log_value": pd.Series(["", pd.NA], dtype="string"),
            "pchembl_value": pd.Series([5.5, pd.NA], dtype="Float64"),
            "standard_value": pd.Series([pd.NA, 50.0], dtype="Float64"),
            "standard_units": pd.Series(["nM", "nM"], dtype="string"),
        }
    )

    augmented, filled = _augment_activity_frame(frame)

    assert list(augmented["salt_chembl_id"]) == ["CHEMBL1", "CHEMBL2"]
    assert str(augmented["salt_chembl_id"].dtype) == "string"
    expected = pd.Series(
        [5.5, float(-math.log10(50e-9))], index=frame.index, dtype="Float64"
    )
    pd.testing.assert_series_equal(
        augmented["log_value"],
        expected,
        check_names=False,
        check_exact=False,
        rtol=1e-6,
    )
    assert {"salt_chembl_id", "log_value"} <= filled


def test_augment_activity_frame__fills_missing_boolean_defaults() -> None:
    frame = pd.DataFrame(
        {
            "approx_cited_activity": pd.Series([pd.NA, True], dtype="boolean"),
            "shuffled_cit": pd.Series([pd.NA, pd.NA], dtype="boolean"),
            "exact_cited_activity": pd.Series([pd.NA, False], dtype="boolean"),
            "higly_correlated_cit": pd.Series([pd.NA, pd.NA], dtype="boolean"),
            "review_doc": pd.Series([pd.NA, pd.NA], dtype="boolean"),
            "rounded_data_citation": pd.Series([pd.NA, pd.NA], dtype="boolean"),
        }
    )

    augmented, filled = _augment_activity_frame(frame)

    expectations = {
        "approx_cited_activity": [False, True],
        "shuffled_cit": [False, False],
        "exact_cited_activity": [False, False],
        "higly_correlated_cit": [False, False],
        "review_doc": [False, False],
        "rounded_data_citation": [False, False],
    }

    for column, values in expectations.items():
        assert str(augmented[column].dtype) == "boolean"
        assert augmented[column].tolist() == values
        assert column in filled


def test_rename_columns__propagates_backfilled_pA_value() -> None:
    frame = pd.DataFrame(
        {
            "log_value": pd.Series([pd.NA, ""], dtype="string"),
            "pchembl_value": pd.Series([6.0, 5.0], dtype="Float64"),
        }
    )

    augmented, _ = _augment_activity_frame(frame)
    renamed = _rename_columns(augmented)

    assert "pA_value" in renamed.columns
    pd.testing.assert_series_equal(
        renamed["pA_value"],
        pd.Series([6.0, 5.0], index=frame.index, dtype="Float64"),
        check_names=False,
    )


def test_load_citation_fraction__missing_file(tmp_path: Path) -> None:
    dictionary_root = tmp_path / "dictionary"
    dictionary_root.mkdir()

    with pytest.raises(ActivityExtendedError, match="citation_fraction.csv not found"):
        _load_citation_fraction(dictionary_root)


def test_load_citation_fraction__typed_columns(activity_resources: Path) -> None:
    dictionary_root = activity_resources / "dictionary"

    frame = _load_citation_fraction(dictionary_root)

    assert list(frame.columns) == [
        "N",
        "K_min_significant",
        "test_used_at_threshold",
        "p_value_at_threshold",
    ]
    assert str(frame["N"].dtype) == "Int64"
    assert str(frame["K_min_significant"].dtype) == "Int64"


def test_load_target_metadata__string_columns(activity_resources: Path) -> None:
    dictionary_root = activity_resources / "dictionary"
    path = dictionary_root / "_target" / "targets_type.csv"

    frame = _load_target_metadata(path)

    assert set(frame.columns) == {
        "target_chembl_id",
        "sortorder.target",
        "multifunctional_enzyme",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "genus",
        "superkingdom",
        "phylum",
        "taxon_id",
        "gene_index",
        "taxon_index",
        "unicellular_organism",
    }
    expected_logical = {"multifunctional_enzyme", "unicellular_organism"}
    expected_integer = {"taxon_id"}
    for column, dtype in frame.dtypes.items():
        if column in expected_logical:
            assert dtype == "boolean"
        elif column in expected_integer:
            assert dtype == "Int64"
        else:
            assert dtype == "string[python]"


def test_load_target_metadata__cp1252_roundtrip(tmp_path: Path) -> None:
    target_path = tmp_path / "targets_type.csv"
    rows = [
        (
            "target_chembl_id,target_sort_order,multifunctional_enzyme,IUPHAR_class,"
            "IUPHAR_subclass,genus,superkingdom,phylum,taxon_id,gene_index,taxon_index,"
            "unicellular_organism\n"
        ),
        (
            "TAR-CP1252,010,TRUE,Clássé,Subclássé,Escherichia,Protéobactéria,Firmicutes,123,GénE,Táxon,FALSE\n"
        ),
    ]
    target_path.write_text("".join(rows), encoding="cp1252")

    frame = _load_target_metadata(target_path)

    assert frame.loc[0, "IUPHAR_class"] == "Clássé"
    assert frame.loc[0, "sortorder.target"] == "010"
    assert frame.loc[0, "unicellular_organism"].item() is False
    assert frame.loc[0, "taxon_id"].item() == 123
    assert frame.loc[0, "multifunctional_enzyme"].item() is True
    assert frame["sortorder.target"].dtype == "string[python]"
    assert frame["unicellular_organism"].dtype == "boolean"
    assert frame["taxon_id"].dtype == "Int64"
    assert frame["multifunctional_enzyme"].dtype == "boolean"


def test_load_target_metadata__derives_unicellular_flag(tmp_path: Path) -> None:
    target_path = tmp_path / "targets_type.csv"
    target_path.write_text(
        (
            "target_chembl_id,target_sort_order,multifunctional_enzyme,IUPHAR_class,"
            "IUPHAR_subclass,genus,superkingdom,phylum,taxon_id,gene_index,taxon_index,"
            "organism_type\n"
            "CHEMBL1,001,FALSE,ClassA,SubclassA,Homo,Eukaryota,Chordata,9606,G1,T1,Multicellular organism\n"
            "CHEMBL2,002,TRUE,ClassB,SubclassB,Escherichia,Bacteria,Proteobacteria,83333,G2,T2,Unicellular organism\n"
        ),
        encoding="utf-8",
    )

    frame = _load_target_metadata(target_path)

    assert list(frame.columns) == list(_TARGET_COLUMNS)
    assert frame.loc[0, "unicellular_organism"].item() is False
    assert frame.loc[1, "unicellular_organism"].item() is True
    assert frame["unicellular_organism"].dtype == "boolean"


def test_resolve_targets_path__uses_override(tmp_path: Path) -> None:
    dictionary_root = tmp_path / "dictionary"
    dictionary_root.mkdir()
    override = tmp_path / "custom_targets.csv"
    override.write_text(
        "target_chembl_id,target_sort_order\nTAR-X,1\n", encoding="utf-8"
    )

    resolved = _resolve_targets_path(dictionary_root, override)

    assert resolved == override


def test_resolve_targets_path__prefers_known_subdir(tmp_path: Path) -> None:
    dictionary_root = tmp_path / "dictionary"
    targets_dir = dictionary_root / "_target"
    targets_dir.mkdir(parents=True)
    expected = targets_dir / "targets_type.csv"
    expected.write_text(
        "target_chembl_id,target_sort_order\nTAR-X,1\n", encoding="utf-8"
    )

    resolved = _resolve_targets_path(dictionary_root, None)

    assert resolved == expected


def test_resolve_targets_path__raises_when_missing(tmp_path: Path) -> None:
    dictionary_root = tmp_path / "dictionary"
    dictionary_root.mkdir()

    with pytest.raises(ActivityExtendedError, match="targets_type.csv not found"):
        _resolve_targets_path(dictionary_root, None)


def test_latest_activity_export__picks_latest_file(tmp_path: Path) -> None:
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    older = search_dir / "output.activity_20230101.csv"
    newer = search_dir / "output.activity_20240101.csv"
    older.write_text("activity_chembl_id\nACT-1\n", encoding="utf-8")
    newer.write_text("activity_chembl_id\nACT-2\n", encoding="utf-8")

    resolved = _latest_activity_export(search_dir)

    assert resolved == newer


def test_latest_activity_export__handles_hidden_temporary_file(tmp_path: Path) -> None:
    search_dir = tmp_path / "exports"
    search_dir.mkdir()
    visible = search_dir / "output.activity_20230101.csv"
    hidden_tmp = search_dir / ".output.activities_20240101.csv.tmp"
    visible.write_text("activity_chembl_id\nACT-1\n", encoding="utf-8")
    hidden_tmp.write_text("activity_chembl_id\nACT-2\n", encoding="utf-8")

    resolved = _latest_activity_export(search_dir)

    assert resolved == hidden_tmp


def test_derive_output_path__normalises_plural_basename(tmp_path: Path) -> None:
    source = tmp_path / ".output.activities_20240101.csv.tmp"

    derived = _derive_output_path(source)

    assert derived == source.with_name("extended.output.activity_20240101.csv")


@pytest.mark.parametrize(
    "source_name",
    [
        "extended.output.activity_20240101.csv",
        "extended.extended.output.activity_20240101.csv.tmp",
    ],
)
def test_derive_output_path__deduplicates_extended_prefix(
    tmp_path: Path, source_name: str
) -> None:
    source = tmp_path / source_name

    derived = _derive_output_path(source)

    assert derived == source.with_name("extended.output.activity_20240101.csv")


def test_derive_output_path__prefixes_custom_exports(tmp_path: Path) -> None:
    source = tmp_path / "chembl_activities_snapshot.csv"

    derived = _derive_output_path(source)

    assert derived.name == "extended.chembl_activities_snapshot.csv"


def test_transform_activity_frame__parses_activity_properties_flags(
    activity_resources: Path,
) -> None:
    dictionary_root = activity_resources / "dictionary"
    export_path = activity_resources / "exports" / "output.activity_20240101.csv"
    frame = pd.read_csv(export_path)

    transformed = _transform_activity_frame(
        frame,
        dictionary_root=dictionary_root,
        targets_override=None,
    )

    assert list(transformed.columns) == list(_FINAL_COLUMN_ORDER)

    flag_map: dict[str, dict[str, bool]] = {}
    for row in frame.itertuples(index=False):
        payload = json.loads(row.activity_properties)
        flag_map[row.activity_chembl_id] = {
            "exact_data_citation": bool(payload["flags"]["exact_data_citation"]),
            "higly_correlated_assay": bool(payload["flags"]["higly_correlated_assay"]),
            "shuffled_assay": bool(payload["flags"]["shuffled_assay"]),
            "review": bool(payload["flags"]["review"]),
            "rounded_data_citation": bool(payload["flags"]["rounded_data_citation"]),
        }

    for row in transformed.itertuples(index=False):
        flags = flag_map[row.activity_chembl_id]
        assert bool(row.exact_data_citation) == flags["exact_data_citation"]
        assert bool(row.higly_correlated_assay) == flags["higly_correlated_assay"]
        assert bool(row.shuffled_assay) == flags["shuffled_assay"]
        assert bool(row.review) == flags["review"]
        assert bool(row.rounded_data_citation) == flags["rounded_data_citation"]

    assert transformed["completed"].tolist() == [
        "1980-08-15",
        "1980-08-15",
        "2008-01-11",
        pd.NA,
    ]
    assert transformed["assay_with_same_target"].tolist() == [3, 3, 2, pd.NA]
    assert transformed["molecule_chembl_id.1"].tolist() == [
        "MOL-1",
        "MOL-1",
        "MOL-2",
        "MOL-3",
    ]
    assert transformed["standard_inchi_skeleton"].tolist() == [
        "InChI=1",
        "InChI=1",
        "InChI=2",
        pd.NA,
    ]

    dtype_map = {
        column: str(transformed[column].dtype) for column in transformed.columns
    }
    assert dtype_map == EXPECTED_DTYPE_MAP

    assert transformed["allosteric"].dtype == pd.BooleanDtype()
    assert transformed["nam"].dtype == pd.BooleanDtype()
    assert transformed["pam"].dtype == pd.BooleanDtype()
    assert transformed["allosteric"].eq(False).all()
    assert transformed["nam"].eq(False).all()
    assert transformed["pam"].eq(False).all()


def test_transform_activity_frame__fills_missing_columns(
    activity_resources: Path,
) -> None:
    dictionary_root = activity_resources / "dictionary"
    frame = pd.DataFrame(
        {
            "activity_id": ["ACT-001"],
            "assay_chembl_id": ["ASSAY-1"],
            "molecule_chembl_id": ["MOL-1"],
            "target_chembl_id": ["TAR-ALPHA"],
            "document_chembl_id": ["DOC-1"],
            "bao_endpoint": ["BAO:000001"],
            "bao_format": ["FMT-1"],
            "standard_type": ["IC50"],
            "standard_value": [8.0],
            "pchembl_value": [5.0],
            "molecule_pref_name": ["Compound One"],
        }
    )

    transformed = _transform_activity_frame(
        frame,
        dictionary_root=dictionary_root,
        targets_override=None,
    )

    assert list(transformed.columns) == list(_FINAL_COLUMN_ORDER)
    row = transformed.iloc[0]
    assert row.activity_chembl_id == "ACT-001"
    assert row.compound_key == "CMPD-1"
    assert row.compound_name == "Compound One"
    assert row.saltform_id == "MOL-1"
    assert row.completed == "1980-08-15"
    assert row.assay_with_same_target == 3
    assert row["molecule_chembl_id.1"] == "MOL-1"
    assert row.standard_inchi_skeleton == "InChI=1"
    assert bool(row.multmol_assay) is False
    assert bool(row.rounded_data_citation) is False
    assert pd.isna(row.original_activity_approx)
    assert pd.isna(row.original_activity_exact)
    assert float(row.pA_value) == pytest.approx(5.0)


def test_apply_multimol_logic__marks_duplicate_multimol_assay(
    activity_resources: Path,
) -> None:
    export_path = activity_resources / "exports" / "output.activity_20240101.csv"
    frame = pd.read_csv(export_path).head(2)
    prepared = _prepare_unknown_chirality(frame)

    result = _apply_multimol_logic(prepared)

    assert result["multmol_assay"].tolist() == [True, True]


def test_select_and_cast__missing_columns_error() -> None:
    frame = pd.DataFrame({"activity_chembl_id": ["ACT-1"]})

    with pytest.raises(ActivityExtendedError, match="missing expected columns"):
        _select_and_cast(frame)


def test_transform_activity_frame__fills_missing_required_columns(
    activity_resources: Path,
) -> None:
    dictionary_root = activity_resources / "dictionary"
    frame = pd.DataFrame(
        {
            "activity_id": ["ACT-1"],
            "assay_chembl_id": ["ASSAY-1"],
            "document_chembl_id": ["DOC-1"],
            "molecule_chembl_id": ["CHEMBL1"],
            "target_chembl_id": ["TAR-ALPHA"],
            "standard_type": ["IC50"],
            "standard_value": [12.3],
            "bao_format": ["Single protein"],
            "pchembl_value": [6.5],
        }
    )

    transformed = _transform_activity_frame(
        frame,
        dictionary_root=dictionary_root,
        targets_override=None,
    )

    assert list(transformed.columns) == list(_FINAL_COLUMN_ORDER)
    assert transformed.loc[0, "saltform_id"] == "CHEMBL1"
    assert float(transformed.loc[0, "pA_value"]) == pytest.approx(6.5)
    assert transformed.loc[0, "compound_name"] is pd.NA
    assert transformed.loc[0, "completed"] == "1980-08-15"
    assert transformed.loc[0, "assay_with_same_target"] == 3
    assert pd.isna(transformed.loc[0, "molecule_chembl_id.1"])
    assert pd.isna(transformed.loc[0, "standard_inchi_skeleton"])


def test_process_activity_extended__writes_expected_payload(
    activity_resources: Path,
    tmp_path: Path,
) -> None:
    tmp_exports = _copytree(activity_resources / "exports", tmp_path / "exports")
    tmp_dictionary = _copytree(
        activity_resources / "dictionary", tmp_path / "dictionary"
    )

    output_path = process_activity_extended(
        search_dir=tmp_exports,
        dictionary_dir=tmp_dictionary,
    )

    assert output_path.name == "extended.output.activity_20240101.csv"

    source_frame = pd.read_csv(tmp_exports / "output.activity_20240101.csv")
    transformed = _transform_activity_frame(
        source_frame,
        dictionary_root=Path(tmp_dictionary),
        targets_override=None,
    )
    assert list(transformed.columns) == list(_FINAL_COLUMN_ORDER)
    dtype_map = {
        column: str(transformed[column].dtype) for column in transformed.columns
    }
    assert dtype_map == EXPECTED_DTYPE_MAP

    result = pd.read_csv(output_path)
    expected = pd.read_csv(
        activity_resources / "expected.extended.output.activity_20240101.csv"
    )

    pd.testing.assert_frame_equal(result, expected)


def test_process_activity_extended__supports_base_dir_alias(
    activity_resources: Path,
    tmp_path: Path,
) -> None:
    tmp_exports = _copytree(activity_resources / "exports", tmp_path / "exports")
    tmp_dictionary = _copytree(
        activity_resources / "dictionary", tmp_path / "dictionary"
    )

    with pytest.warns(DeprecationWarning, match="base_dir"):
        output_path = process_activity_extended(
            base_dir=tmp_exports,
            dictionary_dir=tmp_dictionary,
        )

    assert output_path.name == "extended.output.activity_20240101.csv"


def test_process_activity_extended__raises_when_search_and_base_dir_provided(
    activity_resources: Path,
    tmp_path: Path,
) -> None:
    tmp_exports = _copytree(activity_resources / "exports", tmp_path / "exports")
    tmp_dictionary = _copytree(
        activity_resources / "dictionary", tmp_path / "dictionary"
    )

    with pytest.raises(TypeError, match="both 'search_dir' and 'base_dir'"):
        process_activity_extended(
            search_dir=tmp_exports,
            base_dir=tmp_exports,
            dictionary_dir=tmp_dictionary,
        )


def test_process_activity_extended__raises_for_missing_dictionary(
    activity_resources: Path,
    tmp_path: Path,
) -> None:
    tmp_exports = _copytree(activity_resources / "exports", tmp_path / "exports")
    missing_dictionary = tmp_path / "dictionary"
    missing_dictionary.mkdir()

    with pytest.raises(ActivityExtendedError, match="molecule_hierarchy.csv not found"):
        process_activity_extended(
            search_dir=tmp_exports,
            dictionary_dir=missing_dictionary,
        )


def test_process_activity_extended__uses_explicit_input_path(
    activity_resources: Path,
    tmp_path: Path,
) -> None:
    tmp_exports = _copytree(activity_resources / "exports", tmp_path / "exports")
    tmp_dictionary = _copytree(
        activity_resources / "dictionary", tmp_path / "dictionary"
    )

    original = tmp_exports / "output.activity_20240101.csv"
    explicit = tmp_exports / "chembl_activities_snapshot.csv"
    original.rename(explicit)

    output_path = process_activity_extended(
        search_dir=tmp_exports,
        input_path=explicit,
        dictionary_dir=tmp_dictionary,
    )

    assert output_path.name == "extended.chembl_activities_snapshot.csv"

    result = pd.read_csv(output_path)
    expected = pd.read_csv(
        activity_resources / "expected.extended.output.activity_20240101.csv"
    )

    pd.testing.assert_frame_equal(result, expected)
