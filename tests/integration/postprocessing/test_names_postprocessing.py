from __future__ import annotations

import importlib.util
import json
import sys
import types
from functools import partial
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

import library

from . import assert_csv_header

REPO_ROOT = Path(__file__).resolve().parents[3]


def _load_names_module() -> types.ModuleType:
    module_name = "library.postprocessing.names"
    module = sys.modules.get(module_name)
    if module is not None:
        return module

    package_name = "library.postprocessing"
    package = sys.modules.get(package_name)
    if package is None:
        package = types.ModuleType(package_name)
        package.__path__ = [str(REPO_ROOT / "library" / "postprocessing")]
        sys.modules[package_name] = package
        library.postprocessing = package

    path = REPO_ROOT / "library" / "postprocessing" / "names.py"
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None:
        raise ImportError("Unable to create spec for names module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module

    source = path.read_text(encoding="utf-8")
    fixed = source.replace(
        '"""Helpers for generating target name exports.',
        "Helpers for generating target name exports.",
        1,
    )
    exec(compile(fixed, str(path), "exec"), module.__dict__)

    package.names = module
    return module


names = _load_names_module()


@pytest.mark.unit
def test_is_hidrate__detects_known_annotations() -> None:
    assert names.is_hidrate("Calcium dihydrate") is True
    assert names.is_hidrate("anhydrous sodium chloride") is False


@pytest.mark.unit
def test_remove_hidrate__removes_tokens_and_collapses_spaces() -> None:
    text = "  Magnesium sesquihydrate   "
    assert names.remove_hidrate(text) == "Magnesium"


@pytest.mark.unit
def test_sort_my_list__deduplicates_and_preserves_stable_order() -> None:
    values = ["Beta", "alpha", "Alpha", None, "beta", "Gamma", "Beta"]
    result = names.sort_my_list(values)
    assert result == "alpha|Alpha|Beta|beta|Gamma"


@pytest.mark.unit
def test_reference_smiles__override_priority_over_table(tmp_path: Path) -> None:
    table_path = tmp_path / "Table6.csv"
    table_path.write_text(
        "molecule_chembl_id,canonical_smiles\nCHEMBL1,C1=CC=CC=C1\n",
        encoding="utf-8",
    )

    result = names.reference_SMILES(
        "CHEMBL1",
        overrides={"CHEMBL1": "C[C@H](O)C"},
        reference_path=table_path,
    )

    assert result == "C[C@H](O)C"


@pytest.mark.unit
def test_reference_smiles__missing_table_raises_error(tmp_path: Path) -> None:
    missing = tmp_path / "Table6.csv"
    with pytest.raises(
        names.TargetNamesError, match="Reference SMILES table not found"
    ):
        names.reference_SMILES("CHEMBL1", reference_path=missing)


@pytest.mark.unit
def test_reference_smiles__missing_required_columns_raises(tmp_path: Path) -> None:
    bad_table = tmp_path / "Table6.csv"
    bad_table.write_text("molecule_chembl_id\nCHEMBL1\n", encoding="utf-8")

    with pytest.raises(names.TargetNamesError, match="missing required columns"):
        names.reference_SMILES("CHEMBL1", reference_path=bad_table)


@pytest.mark.unit
def test_component_rows__hydrates_salts_and_structures(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference_path = tmp_path / "Table6.csv"
    reference_path.write_text(
        "molecule_chembl_id,canonical_smiles\nCHEMBL100,C1CCO1\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(names, "_REFERENCE_CACHE", {})
    monkeypatch.setattr(
        names,
        "reference_SMILES",
        partial(names.reference_SMILES, reference_path=reference_path),
    )

    components = [
        {
            "molecule_chembl_id": "CHEMBL100",
            "component_name": "Alpha hemihydrate",
            "component_synonyms": [
                {"synonym": "Alpha hemihydrate"},
                {"synonym": "Alpha"},
            ],
            "molecule_structures": {
                "standard_inchi": "InChI=1",
                "standard_inchi_stereo": "InChI=1S",
            },
            "molecule_properties": {"mw_freebase": 123.4},
            "unknown_chirality": True,
            "contrion": ["B", "A"],
            "component_type": "Protein",
        },
        {
            "component_chembl_id": "CHEMBL200",
            "component_name": "Beta hydrate",
            "component_synonyms": [{"synonym": "Beta"}],
            "molecule_structures": {
                "canonical_smiles": "C",
                "standard_inchi": "InChI=2",
                "standard_inchi_stereo": "InChI=2S",
            },
            "molecule": {
                "molecule_structures": {
                    "canonical_smiles": "C",
                    "standard_inchi": "InChI=2",
                    "standard_inchi_stereo": "InChI=2S",
                },
                "molecule_hierarchy": {"parent_molecule_chembl_id": "CHEMBL999"},
            },
            "molecule_properties": {"mw_freebase": 456.7},
            "unknown_chirality": False,
            "contrion": "C|D",
            "active_component": "",  # forces fallback to name
            "component_type": "Complex",
        },
    ]

    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL42",
                "target_components": json.dumps(components),
            }
        ]
    )

    rows = names._component_rows(frame)

    assert len(rows) == 2

    first, second = rows

    assert first["canonical_smiles"] == "C1CCO1"
    assert first["is_hydrate"] == "Y"
    assert first["is_saltform"] == "N"
    assert first["unknown_chirality"] == "Y"
    assert first["standard_inchi_skeleton"] == "InChI=1"
    assert first["standard_inchi_stereo"] == "InChI=1S"
    assert first["contrion"] == "A|B"
    assert first["active_component"] == "Alpha"
    assert first["active_component_type"] == "Protein"

    assert second["canonical_smiles"] == "C"
    assert second["salt_chembl_id"] == "CHEMBL999"
    assert second["is_saltform"] == "Y"
    assert second["is_hydrate"] == "Y"
    assert second["active_component"] == "Beta"
    assert second["active_component_type"] == "Complex"


@pytest.mark.unit
def test_component_rows__missing_identifier_raises_error() -> None:
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL123",
                "target_components": json.dumps([{"component_synonyms": []}]),
            }
        ]
    )

    with pytest.raises(names.TargetNamesError, match="missing molecule_chembl_id"):
        names._component_rows(frame)


@pytest.mark.unit
def test_process_target_names_helper__writes_byte_identical_output(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "output.target_20250101.csv"
    data = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL100",
                "uniprot_id_primary": "P11111",
                "pref_name": "Alpha",
                "protein_name_canonical": "Alpha canonical",
                "synonyms": "Alpha|Beta|alpha",
                "protein_synonym_list": "",
                "gtop_synonyms": "Alpha GPCR",
                "gene_symbol": "ALPHA",
                "gene_symbol_list": "ALPHA|Alpha",
                "isoform_names": "Isoform A|Isoform B",
                "isoform_synonyms": "",
                "recommendedName": "Alpha recommended",
                "secondaryAccessionNames": "SEC1|SEC2",
                "contrion": "X|Y",
                "active_component_type": "protein",
            },
            {
                "target_chembl_id": "CHEMBL200",
                "uniprot_id_primary": "Q22222",
                "pref_name": "",
                "protein_name_canonical": "",
                "synonyms": "-|-",
                "protein_synonym_list": "Gamma",
                "gtop_synonyms": "",
                "gene_symbol": "",
                "gene_symbol_list": "",
                "isoform_names": "",
                "isoform_synonyms": "Gamma iso",
                "recommendedName": "",
                "secondaryAccessionNames": "",
                "contrion": "",
                "active_component_type": "complex",
            },
        ]
    )
    data.to_csv(input_path, index=False)

    result = names.process_target_names(input_path)
    output_path = Path(result["path"])

    produced = pd.read_csv(output_path, dtype="string")
    assert_csv_header(output_path, names.TARGET_NAMES_COLUMNS)

    expected_records = [
        ("CHEMBL100", "P11111", "Alpha", "chembl_preferred", "pref_name"),
        (
            "CHEMBL100",
            "P11111",
            "Alpha canonical",
            "uniprot_canonical",
            "protein_name_canonical",
        ),
        ("CHEMBL100", "P11111", "Alpha GPCR", "gtop_synonym", "gtop_synonyms"),
        (
            "CHEMBL100",
            "P11111",
            "Alpha recommended",
            "uniprot_recommended",
            "recommendedName",
        ),
        ("CHEMBL100", "P11111", "Alpha", "chembl_synonym", "synonyms"),
        ("CHEMBL100", "P11111", "Beta", "chembl_synonym", "synonyms"),
        ("CHEMBL100", "P11111", "alpha", "chembl_synonym", "synonyms"),
        ("CHEMBL100", "P11111", "ALPHA", "gene_symbol", "gene_symbol"),
        ("CHEMBL100", "P11111", "ALPHA", "gene_symbol", "gene_symbol_list"),
        ("CHEMBL100", "P11111", "Alpha", "gene_symbol", "gene_symbol_list"),
        ("CHEMBL100", "P11111", "Isoform A", "isoform_name", "isoform_names"),
        ("CHEMBL100", "P11111", "Isoform B", "isoform_name", "isoform_names"),
        ("CHEMBL100", "P11111", "SEC1", "uniprot_secondary", "secondaryAccessionNames"),
        ("CHEMBL100", "P11111", "SEC2", "uniprot_secondary", "secondaryAccessionNames"),
        ("CHEMBL200", "Q22222", "Gamma", "uniprot_synonym", "protein_synonym_list"),
        ("CHEMBL200", "Q22222", "Gamma iso", "isoform_synonym", "isoform_synonyms"),
    ]

    expected = pd.DataFrame(
        expected_records,
        columns=names.TARGET_NAMES_COLUMNS,
    ).astype("string")

    expected = expected.sort_values(
        by=["target_chembl_id", "name", "name_type"], kind="mergesort"
    ).reset_index(drop=True)

    pdt.assert_frame_equal(produced, expected)

    expected_path = tmp_path / "expected.csv"
    names.write_csv_deterministic(
        expected,
        expected_path,
        col_order=names.TARGET_NAMES_COLUMNS,
        key_cols=["target_chembl_id"],
        encoding="utf-8",
        sep=",",
        cfg=None,
    )
    assert output_path.read_bytes() == expected_path.read_bytes()


@pytest.mark.unit
def test_process_target_names_helper__summary_counts(tmp_path: Path) -> None:
    input_path = tmp_path / "output.target_20250101.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL100",
                "uniprot_id_primary": "P11111",
                "pref_name": "Alpha",
                "contrion": "X|Y",
                "active_component_type": "protein",
            },
            {
                "target_chembl_id": "CHEMBL200",
                "uniprot_id_primary": "Q22222",
                "pref_name": "Beta",
                "contrion": "",
                "active_component_type": "complex",
            },
        ]
    )
    frame.to_csv(input_path, index=False)

    result = names.process_target_names(input_path)
    assert_csv_header(Path(result["path"]), names.TARGET_NAMES_COLUMNS)

    summary = result["summary"]
    assert summary["rows_before"] == 2
    assert summary["rows_after"] == 2
    assert summary["contrion_unique"] == 2
    assert summary["contrion_non_empty"] == 1
    assert summary["contrion_total"] == 2
    assert summary["active_component_type"] == {"protein": 1, "complex": 1}


@pytest.mark.unit
def test_process_target_names_helper__verbose_logs_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    input_path = tmp_path / "output.target_20250101.csv"
    pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL100",
                "uniprot_id_primary": "P11111",
                "pref_name": "Alpha",
            }
        ]
    ).to_csv(input_path, index=False)

    events: list[tuple[str, dict[str, object]]] = []

    def _record(event: str, **extra: object) -> None:
        events.append((event, extra))

    monkeypatch.setattr(names.logger, "info", _record)

    result = names.process_target_names(input_path, verbose=True)
    assert events, "Expected logger.info to be called when verbose=True"
    event, payload = events[0]
    assert event == "target_names_helper_done"
    assert payload["path"] == str(Path(result["path"]))
    assert payload["rows"] == 1
    assert Path(result["path"]).name.startswith("names.output.target_20250101.csv")
    assert_csv_header(Path(result["path"]), names.TARGET_NAMES_COLUMNS)


@pytest.mark.unit
def test_process_target_names_helper__stable_sorting_of_duplicates(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "output.target_20250101.csv"
    pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL300",
                "uniprot_id_primary": "P33333",
                "gene_symbol": "BRCA1",
                "gene_symbol_list": "BRCA1|BRCA1-alt",
            }
        ]
    ).to_csv(input_path, index=False)

    result = names.process_target_names(input_path)
    output_path = Path(result["path"])
    produced = pd.read_csv(output_path, dtype="string")
    assert_csv_header(output_path, names.TARGET_NAMES_COLUMNS)

    subset = produced[produced["name"] == "BRCA1"]
    assert list(subset["source_column"]) == ["gene_symbol", "gene_symbol_list"]


def test_process_target_names__custom_output_dir(tmp_path: Path) -> None:
    output_dir = tmp_path / "custom-output"
    output_dir.mkdir()

    older = output_dir / "output.target_20240101.csv"
    newer = output_dir / "output.target_20250303.csv"

    pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_OLD",
                "uniprot_id_primary": "POLD",
                "pref_name": "Legacy",
            }
        ]
    ).to_csv(older, index=False)

    pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_NEW",
                "uniprot_id_primary": "PNEW",
                "pref_name": "Latest",
                "contrion": "A|B",
                "active_component_type": "protein",
            }
        ]
    ).to_csv(newer, index=False)

    result = names.process_target_names(output_dir=output_dir)
    output_path = Path(result["path"])

    assert output_path.parent == output_dir
    assert output_path.name == "names.output.target_20250303.csv"

    produced = pd.read_csv(output_path, dtype="string")
    assert_csv_header(output_path, names.TARGET_NAMES_COLUMNS)

    assert (produced["target_chembl_id"] == "CHEMBL_NEW").all()
    assert "Latest" in produced["name"].tolist()

    summary = result["summary"]
    assert summary["rows_before"] == 1
    assert summary["contrion_total"] == 2


@pytest.mark.unit
def test_ensure_columns__missing_required_columns_raise() -> None:
    frame = pd.DataFrame({"target_components": ["[]"]})
    with pytest.raises(names.TargetNamesError, match="missing required columns"):
        names._ensure_columns(frame, ["target_chembl_id", "target_components"])
