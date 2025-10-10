from __future__ import annotations

import csv
from pathlib import Path

from library.postprocessing.activity_extended import _load_target_metadata


def _write_target_csv(path: Path, rows: list[dict[str, object]]) -> None:
    columns = [
        "target_chembl_id",
        "target_sort_order",
        "multifunctional_enzyme",
        "IUPHAR_class",
        "IUPHAR_subclass",
        "genus",
        "superkingdom",
        "phylum",
        "taxon_id",
        "gene_index",
        "taxon_index",
    ]
    if rows and "Organism Type" in rows[0]:
        columns.append("Organism Type")
    elif rows and "Unicellular organism" in rows[0]:
        columns.append("Unicellular organism")
    else:
        columns.append("unicellular_organism")

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _base_row(**overrides: object) -> dict[str, object]:
    base = {
        "target_chembl_id": "CHEMBL1",
        "target_sort_order": "001",
        "multifunctional_enzyme": "FALSE",
        "IUPHAR_class": "Class",
        "IUPHAR_subclass": "Subclass",
        "genus": "Genus",
        "superkingdom": "Bacteria",
        "phylum": "Firmicutes",
        "taxon_id": 9606,
        "gene_index": "GENE",
        "taxon_index": "9606",
    }
    base.update(overrides)
    return base


def test_load_target_metadata__derives_unicellular_from_type_alias(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "targets.csv"
    _write_target_csv(
        csv_path,
        [
            _base_row(**{"Organism Type": "Unicellular organism"}),
            _base_row(
                target_chembl_id="CHEMBL2",
                **{"Organism Type": "Multicellular organism"},
            ),
        ],
    )

    frame = _load_target_metadata(csv_path)

    assert list(frame.columns) == [
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
    ]
    assert frame["unicellular_organism"].tolist() == [True, False]


def test_load_target_metadata__normalises_unicellular_column_name(
    tmp_path: Path,
) -> None:
    csv_path = tmp_path / "targets.csv"
    _write_target_csv(
        csv_path,
        [
            _base_row(**{"Unicellular organism": True}),
            _base_row(target_chembl_id="CHEMBL2", **{"Unicellular organism": False}),
        ],
    )

    frame = _load_target_metadata(csv_path)

    assert "unicellular_organism" in frame.columns
    assert str(frame["unicellular_organism"].dtype) in {"boolean", "bool"}
    assert frame["unicellular_organism"].tolist() == [True, False]
