from pathlib import Path

import pandas as pd

from library.postprocessing import names


def test_process_target_names__creates_name_output_file(tmp_path):
    input_path = tmp_path / "output.target_20240101.csv"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_TARGET",
                "uniprot_id_primary": "P12345",
                "pref_name": "Alpha",
                "synonyms": "Beta|Gamma",
                "gtop_synonyms": "Delta",
                "gene_symbol": "GENE1",
                "gene_symbol_list": "GENE1|GENE2",
                "isoform_names": "Isoform A",
                "isoform_synonyms": "IsoSyn|IsoSyn2",
                "secondaryAccessionNames": "Secondary",
                "contrion": "Token1|Token2",
                "active_component_type": "Protein",
            }
        ],
        dtype="string",
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    result = names.process_target_names(input_path, verbose=False)

    assert isinstance(result, dict)
    output_path = Path(result["path"])
    assert output_path.name == "names.output.target_20240101.csv"
    assert output_path.exists()

    exported = pd.read_csv(output_path, dtype=str)
    assert not exported.empty
    assert (exported["target_chembl_id"] == "CHEMBL_TARGET").all()
    assert "Alpha" in exported["name"].values


def test_process_target_names__normalises_tmp_suffix(tmp_path: Path) -> None:
    input_path = tmp_path / "output.target_20250101.csv.tmp"
    frame = pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_TARGET",
                "uniprot_id_primary": "P54321",
                "pref_name": "Beta",
                "gene_symbol": "GENE2",
            }
        ],
        dtype="string",
    )
    frame.to_csv(input_path, index=False, encoding="utf-8")

    result = names.process_target_names(input_path, verbose=False)

    output_path = Path(result["path"])
    assert output_path.name == "names.output.target_20250101.csv"
    assert output_path.exists()


def test_process_target_names__writes_to_custom_output_dir(tmp_path: Path) -> None:
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    input_path = input_dir / "output.target_20260101.csv"
    pd.DataFrame(
        [
            {
                "target_chembl_id": "CHEMBL_CUSTOM",
                "uniprot_id_primary": "PCUSTOM",
                "pref_name": "Custom",
            }
        ]
    ).to_csv(input_path, index=False)

    output_dir = tmp_path / "names"
    result = names.process_target_names(input_path, output_dir=output_dir, verbose=False)

    output_path = Path(result["path"])
    assert output_path.parent == output_dir
    assert output_path.name == "names.output.target_20260101.csv"
    assert output_path.exists()
