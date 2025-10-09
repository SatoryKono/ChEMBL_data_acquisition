from __future__ import annotations

import csv
import importlib.util
import sys
import types
import warnings
from pathlib import Path


import importlib.util
import sys
import types
from typing import Any

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

import library
from library.config import Config
from library.pipelines.target import cellularity, helpers, multifunctional
from library.postprocessing import names as target_names
from library.postprocessing.names import process_target_names
from library.postprocessing.target import isoform
 

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_postprocessing_cellularity_module():
    module = sys.modules.get("library.postprocessing.target.cellularity")
    if module is not None:
        return module

    post_pkg = sys.modules.get("library.postprocessing")
    if post_pkg is None:
        post_pkg = types.ModuleType("library.postprocessing")
        post_pkg.__path__ = [str(REPO_ROOT / "library" / "postprocessing")]
        sys.modules["library.postprocessing"] = post_pkg
        setattr(library, "postprocessing", post_pkg)

    target_pkg = sys.modules.get("library.postprocessing.target")
    if target_pkg is None:
        target_pkg = types.ModuleType("library.postprocessing.target")
        target_pkg.__path__ = [
            str(REPO_ROOT / "library" / "postprocessing" / "target")
        ]
        sys.modules["library.postprocessing.target"] = target_pkg
    if not hasattr(sys.modules["library.postprocessing"], "target"):
        setattr(sys.modules["library.postprocessing"], "target", target_pkg)

    spec = importlib.util.spec_from_file_location(
        "library.postprocessing.target.cellularity",
        REPO_ROOT / "library" / "postprocessing" / "target" / "cellularity.py",
    )
    if spec is None or spec.loader is None:
        raise ImportError("Unable to load postprocessing cellularity module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


postprocessing_cellularity = _load_postprocessing_cellularity_module()
from library.schemas.targets import CELLULARITY_COLUMN_NAME, TARGETS_COLUMN_ORDER

INPUT_FILE = "target_postprocess_power_query_input.csv"
EXPECTED_FILE = "target_postprocess_power_query_expected.csv"

_ROOT = Path(__file__).resolve().parents[2]

if "library.postprocessing" not in sys.modules:
    postprocessing_pkg = types.ModuleType("library.postprocessing")
    postprocessing_pkg.__path__ = [str(_ROOT / "library/postprocessing")]
    sys.modules["library.postprocessing"] = postprocessing_pkg

if "library.postprocessing.target" not in sys.modules:
    target_pkg = types.ModuleType("library.postprocessing.target")
    target_pkg.__path__ = [str(_ROOT / "library/postprocessing/target")]
    sys.modules["library.postprocessing.target"] = target_pkg

_SPEC = importlib.util.spec_from_file_location(
    "library.postprocessing.target.cellularity",
    _ROOT / "library/postprocessing/target/cellularity.py",
)
assert _SPEC is not None and _SPEC.loader is not None
_cellularity_module = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _cellularity_module
_SPEC.loader.exec_module(_cellularity_module)
Cellularity = _cellularity_module.Cellularity


@pytest.mark.unit
def test_read_snapshot__validates_required_columns(snapshot_resource: Path) -> None:
    path = snapshot_resource / INPUT_FILE
    with pytest.raises(ValueError, match="missing required columns"):
        helpers.read_snapshot(path, columns=["target_chembl_id", "missing"])


@pytest.mark.unit
def test_postprocess_target_table__matches_power_query_snapshot(
    snapshot_resource: Path,
) -> None:
    input_path = snapshot_resource / INPUT_FILE
    expected_path = snapshot_resource / EXPECTED_FILE

    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)
    expected = pd.read_csv(expected_path, dtype=str, keep_default_na=False).astype(
        "string"
    )

    result = helpers.postprocess_target_table(frame)

    assert list(result.columns) == TARGETS_COLUMN_ORDER
    pdt.assert_frame_equal(result.reset_index(drop=True), expected)


@pytest.mark.unit
def test_postprocess_target_table_file__writes_expected_bytes(
    tmp_path: Path, cfg: Config, snapshot_resource: Path
) -> None:
    input_path = snapshot_resource / INPUT_FILE
    expected_path = snapshot_resource / EXPECTED_FILE

    working_input = tmp_path / INPUT_FILE
    working_output = tmp_path / "targets_final.csv"
    working_input.write_bytes(input_path.read_bytes())

    output_path = helpers.postprocess_target_table_file(
        working_input, working_output, cfg=cfg
    )

    assert output_path == working_output
    # Normalise EOL across platforms
    assert working_output.read_bytes().replace(b"\r\n", b"\n") == expected_path.read_bytes().replace(b"\r\n", b"\n")

    result = pd.read_csv(working_output, dtype=str, keep_default_na=False).astype(
        "string"
    )
    expected = pd.read_csv(expected_path, dtype=str, keep_default_na=False).astype(
        "string"
    )
    pdt.assert_frame_equal(result, expected)


@pytest.mark.unit
def test_process_target_names__strips_bom_headers(tmp_path: Path) -> None:
    """Ensure UTF-8 BOM headers are stripped so downstream helpers continue to work."""

    # UTF-8 with BOM exports would previously keep the BOM in ``target_chembl_id``
    # which meant downstream helpers could not find required columns.
    csv_path = tmp_path / "output.target_20240101.csv"
    with csv_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "target_chembl_id",
                "pref_name",
                "synonyms",
                "contrion",
                "active_component_type",
            ]
        )
        writer.writerow([
            "CHEMBL_TGT",
            "Example Compound",
            "Example Alias|Example Compound",
            "Na|Cl",
            "protein",
        ])

    result_info = process_target_names(csv_path)
    output_path = Path(result_info["path"])

    assert output_path.name == "names.output.target_20240101.csv"
    assert output_path.exists()
    assert result_info["summary"]["rows_after"] == 3

    result = pd.read_csv(output_path, dtype=str, keep_default_na=False).astype("string")
    assert list(result.columns) == list(target_names.TARGET_NAMES_COLUMNS)
    assert set(result["target_chembl_id"]) == {"CHEMBL_TGT"}


@pytest.mark.unit
def test_cellularity_classification__matches_examples() -> None:
    records = [
        {
            "genus": "Homo",
            "lineage_superkingdom": "Eukaryota",
            "lineage_phylum": "Chordata",
            "lineage_class": "Mammalia",
            "expected": cellularity.DEFAULT_RULES.label_multicellular(),
        },
        {
            "genus": "Candida",
            "lineage_superkingdom": "Eukaryota",
            "lineage_phylum": "Ascomycota",
            "lineage_class": "Saccharomycetes",
            "expected": cellularity.DEFAULT_RULES.label_unicellular(),
        },
        {
            "genus": "Influenza A virus",
            "lineage_superkingdom": "Viruses",
            "lineage_phylum": "Negarnaviricota",
            "lineage_class": "Insthoviricetes",
            "expected": cellularity.DEFAULT_RULES.label_viral(),
        },
    ]

    for record in records:
        label = cellularity.classify_record(record)
        assert label == record["expected"]


@pytest.mark.unit
def test_add_cellularity_smart__fills_column(snapshot_resource: Path) -> None:
    input_path = snapshot_resource / INPUT_FILE
    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)

    enriched = cellularity.add_cellularity_smart(
        frame, output_col=CELLULARITY_COLUMN_NAME
    )

    assert CELLULARITY_COLUMN_NAME in enriched.columns
    assert list(enriched[CELLULARITY_COLUMN_NAME]) == [
        cellularity.DEFAULT_RULES.label_multicellular(),
        cellularity.DEFAULT_RULES.label_unicellular(),
        cellularity.DEFAULT_RULES.label_viral(),
    ]


@pytest.mark.unit
def test_add_cellularity_smart__missing_tax_ids_skip_fetch() -> None:
    calls: list[tuple[object, object | None]] = []

    def _stub_fetcher(tax_id: object, email: str | None) -> list[str]:
        calls.append((tax_id, email))
        return ["unexpected"]

    frame = pd.DataFrame(
        {
            "tax_id": pd.Series([pd.NA, np.nan], dtype="object"),
            "lineage_superkingdom": [pd.NA, pd.NA],
            "lineage_phylum": [pd.NA, pd.NA],
        }
    )

    enriched = postprocessing_cellularity.add_cellularity_smart(
        frame,
        tax_id_column="tax_id",
        superkingdom_column="lineage_superkingdom",
        phylum_column="lineage_phylum",
        fetcher=_stub_fetcher,
    )

    assert calls == []
    assert enriched["cellularity"].tolist() == ["ambiguous", "ambiguous"]


@pytest.mark.unit
def test_append_multifunctional_flag__detects_keywords(snapshot_resource: Path) -> None:
    input_path = snapshot_resource / INPUT_FILE
    frame = pd.read_csv(input_path, dtype=str, keep_default_na=False)

    result = multifunctional.append_multifunctional_flag(frame)

    assert "multifunctional_enzyme" in result.columns
    assert result["multifunctional_enzyme"].dtype == "boolean"
    assert result["multifunctional_enzyme"].tolist() == [True, False, False]


@pytest.mark.unit
def test_classify_by_fetch__trailing_whitespace_lineage_remains_ambiguous() -> None:
    lineage = [
        "Viruses ",
        "Duplodnaviria ",
        "Herpesvirales ",
        "Herpesviridae ",
        "Varicellovirus ",
        "Human alphaherpesvirus 3 ",
    ]

    def _mock_fetcher(tax_id: object, email: str | None) -> list[str]:
        assert tax_id == "10335"
        return list(lineage)

    classifier = Cellularity(fetcher=_mock_fetcher)

    names = classifier.get_lineage_names("10335")
    assert names[-1].endswith(" ")

    result = classifier.classify_by_fetch("10335")

    assert result == "ambiguous"
def test_classify_by_fetch__trailing_spaces_remain_ambiguous() -> None:
    lineage = [
        "Cellular organisms",
        "Eukaryota ",
        "Chordata ",
        "Homo sapiens   ",
    ]

    def _fetcher(tax_id: Any, email: str | None) -> list[str]:
        return list(lineage)

    classifier = Cellularity(fetcher=_fetcher)

    result = classifier.classify_by_fetch("9606")

    assert result == "ambiguous"
    assert classifier.get_lineage_names("9606") == lineage


@pytest.mark.unit
def test_get_lineage_names__numeric_tax_id_returns_lineage() -> None:
    lineage = ["Cellular organisms", "Eukaryota", "Chordata", "Homo sapiens"]
    calls: list[tuple[Any, str | None]] = []

    def _fetcher(tax_id: Any, email: str | None) -> list[str]:
        calls.append((tax_id, email))
        return list(lineage)

    classifier = Cellularity(fetcher=_fetcher)

    result = classifier.get_lineage_names(9606)

    assert result == lineage
    assert calls == [(9606, None)]


@pytest.mark.unit
def test_isoform_split_pipes__trims_and_discards_empty() -> None:
    assert isoform._split_pipes("") == []
    assert isoform._split_pipes("a|b") == ["a", "b"]
    assert isoform._split_pipes(" a | | b ") == ["a", "b"]


@pytest.mark.unit
def test_isoform_make_triples__pads_shorter_lists() -> None:
    triples = isoform._make_triples(
        ["name1", "name2"],
        ["id1"],
        ["syn1", "syn2", "syn3"],
    )

    assert triples == [
        {"name": "name1", "id": "id1", "synonym": "syn1"},
        {"name": "name2", "id": None, "synonym": "syn2"},
        {"name": None, "id": None, "synonym": "syn3"},
    ]


@pytest.mark.unit
def test_isoform_normalisation__only_names_and_synonyms_lowercased() -> None:
    frame = pd.DataFrame(
        {
            "isoform_synonyms": ["Alpha|Beta"],
            "isoform_names": ["Alpha|Beta"],
            "isoform_ids": ["EnSp0001|ensp0002"],
            "uniprot_id_primary": ["Q10000"],
            "target_chembl_id": ["CHEMBL100"],
        }
    )

    transform = isoform._transform(frame)

    assert all(name == name.lower() for name in transform.result["name"])
    assert "EnSp0001" in transform.result["id"].tolist()
    assert "ensp0002" in transform.result["id"].tolist()


@pytest.mark.unit
def test_isoform_tokenise_synonym__pde3a_alpha() -> None:
    tokens = isoform._tokenize_synonym("PDE3A:alpha")

    assert tokens == ["pde3a", "3a", "alpha"]


@pytest.mark.unit
def test_isoform_name_filter__drops_empty_and_na(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    assert "" not in transform.result["name"].tolist()
    assert "n/a" not in transform.result["name"].tolist()
    assert "none" not in transform.result["name"].tolist()


@pytest.mark.unit
def test_isoform_synonym_tokens__preserve_identifier(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)
    synonyms = transform.result[transform.result["name"].isin(["pde3a", "3a", "alpha"])]

    assert synonyms[synonyms["name"] == "pde3a"]["id"].tolist() == ["ENSP0002"]
    assert synonyms[synonyms["name"] == "3a"]["id"].tolist() == ["ENSP0002"]
    alpha_ids = synonyms[synonyms["name"] == "alpha"]["id"].tolist()
    assert "ENSP0002" in alpha_ids


@pytest.mark.unit
def test_isoform_stage1_dedup__matches_table_distinct(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    expected = transform.combined.drop_duplicates(
        subset=["id", "name", "target_chembl_id", "uniprot_id_primary"],
        keep="first",
    ).reset_index(drop=True)

    pdt.assert_frame_equal(transform.dedup_stage1.reset_index(drop=True), expected)


@pytest.mark.unit
def test_isoform_sorted_stage__uses_mergesort(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    expected = transform.dedup_stage1.sort_values(
        by=["uniprot_id_primary", "id"],
        kind="mergesort",
        na_position="first",
    ).reset_index(drop=True)

    pdt.assert_frame_equal(transform.sorted_stage.reset_index(drop=True), expected)


@pytest.mark.unit
def test_isoform_stage2_dedup__matches_drop_duplicates(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    expected = transform.sorted_stage.drop_duplicates(
        subset=["id", "target_chembl_id", "name"],
        keep="first",
    ).reset_index(drop=True)

    pdt.assert_frame_equal(transform.dedup_stage2.reset_index(drop=True), expected)


@pytest.mark.unit
def test_isoform_final_dedup__matches_last_distinct(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_duplicates.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    expected = transform.dedup_stage2.drop_duplicates(
        subset=["id", "name"],
        keep="first",
    ).reset_index(drop=True)

    pdt.assert_frame_equal(transform.result.reset_index(drop=True), expected)


@pytest.mark.unit
def test_isoform_output_columns__match_spec(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_minimal.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)

    assert list(transform.result.columns) == list(isoform._OUTPUT_COLUMNS)


@pytest.mark.unit
def test_isoform_make_triples_fixture__pads_missing_values(snapshot_resource: Path) -> None:
    frame = pd.read_csv(
        snapshot_resource / "target_isoform_length_mismatch.csv",
        dtype=str,
        keep_default_na=False,
    )

    transform = isoform._transform(frame)
    rows = transform.combined[transform.combined["target_chembl_id"] == "CHEMBL11"]

    name5_row = rows[rows["name"] == "name5"]
    assert not name5_row.empty
    assert name5_row["id"].apply(lambda value: pd.isna(value) or value == "").all()

    syn5_row = rows[rows["name"] == "syn5"]
    assert not syn5_row.empty
    assert syn5_row["id"].apply(lambda value: pd.isna(value) or value == "").all()


@pytest.mark.unit
def test_isoform_process_targets__writes_isoform_prefix(tmp_path: Path, snapshot_resource: Path) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    input_path = tmp_path / "output.target_20250101.csv"
    input_path.write_bytes(source.read_bytes())

    output_path = isoform.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.output.target_20250101.csv"
    assert output_path.parent == input_path.parent


@pytest.mark.unit
def test_isoform_process_targets__strips_csv_normalized_suffix(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    input_path = tmp_path / "output.targets_20250101.csv_normalized"
    input_path.write_bytes(source.read_bytes())

    match = "does not use the canonical '.csv' extension"
    with pytest.warns(UserWarning, match=match):
        output_path = isoform.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.output.targets_20250101.csv"
    assert output_path.parent == input_path.parent


@pytest.mark.unit
def test_isoform_process_targets__accepts_explicit_custom_name(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    input_path = tmp_path / "custom_export.csv"
    input_path.write_bytes(source.read_bytes())

    match = "does not match the canonical target export naming conventions"
    with pytest.warns(UserWarning, match=match):
        output_path = isoform.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.custom_export.csv"
    assert output_path.parent == input_path.parent


@pytest.mark.unit
def test_isoform_process_targets__warns_on_non_csv_extension(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    source = snapshot_resource / "target_isoform_minimal.csv"
    input_path = tmp_path / "custom_export.txt"
    input_path.write_bytes(source.read_bytes())

    match = "does not use the canonical '.csv' extension"
    with pytest.warns(UserWarning, match=match):
        output_path = isoform.process_targets(str(input_path), verbose=False)

    assert output_path.name == "isoform.custom_export.txt"
    assert output_path.parent == input_path.parent
