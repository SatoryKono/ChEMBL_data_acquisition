from __future__ import annotations

import csv
from pathlib import Path

from library.postprocessing.target.isoform import process_targets


def _write_isoform_input(path: Path) -> None:
    columns = (
        "isoform_synonyms",
        "isoform_names",
        "isoform_ids",
        "uniprot_id_primary",
        "target_chembl_id",
    )
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        writer.writerow(
            [
                "synonym-a|synonym-b",
                "Isoform Alpha",
                "ID-001",
                "P12345",
                "CHEMBL123",
            ]
        )


def test_process_targets__normalises_temp_export_basename(tmp_path: Path) -> None:
    input_path = tmp_path / ".output.targets_20251005.csv.tmp"
    _write_isoform_input(input_path)

    output_path = process_targets(input_csv=str(input_path), verbose=False)

    expected = tmp_path / "isoform.output.targets_20251005.csv"
    assert output_path == expected
    assert output_path.exists()


def test_process_targets__normalises_normalized_extension(tmp_path: Path) -> None:
    input_path = tmp_path / "output.targets_20251005.csv_normalized"
    _write_isoform_input(input_path)

    output_path = process_targets(input_csv=str(input_path), verbose=False)

    expected = tmp_path / "isoform.output.targets_20251005.csv"
    assert output_path == expected
    assert output_path.exists()
