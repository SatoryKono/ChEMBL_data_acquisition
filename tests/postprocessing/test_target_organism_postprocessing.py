from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.postprocessing.target import postprocess_target_table
from library.postprocessing.target.cellularity import Cellularity
from library.postprocessing.target.helpers import normalize_text
from library.postprocessing.target.multifunctional import compute_multifunctional_table


@pytest.mark.unit
def test_normalize_text__mirrors_power_query() -> None:
    assert normalize_text(None) == ""
    assert normalize_text("  Mixed Case ") == "mixed case"
    assert normalize_text(123) == "123"


@pytest.mark.unit
def test_classify_by_lineage__animal_phyly_returns_multicellular() -> None:
    assert Cellularity.classify_by_lineage("Eukaryota", "Chordata") == "multicellular"
    assert Cellularity.classify_by_lineage("Bacteria", "Firmicutes") == "unicellular"
    assert Cellularity.classify_by_lineage("Viruses", None) == "acellular (virus)"


@pytest.mark.unit
def test_add_cellularity_smart__falls_back_to_fetch(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str | None, str | None]] = []

    def fake_fetch(cls, tax_id, email=None):
        calls.append((tax_id, email))
        return "acellular (virus)"

    monkeypatch.setattr(Cellularity, "classify_by_fetch", classmethod(fake_fetch))

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL"],
            "uniprot_id_primary": ["P1"],
            "organism": ["Virus"],
            "taxon_id": ["1234"],
            "lineage_superkingdom": [""],
            "lineage_phylum": [""],
            "lineage_class": [""],
        }
    )

    classified = Cellularity.add_cellularity_smart(
        frame,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
    )
    assert classified.loc[0, "cellularity"] == "acellular (virus)"
    assert calls == [("1234", None)]


@pytest.mark.unit
def test_compute_multifunctional_table__deduplicates_prefixes() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "reaction_ec_numbers": ["2.7.11.1|2.7.11.2|3.1.4.4"],
        }
    )
    result = compute_multifunctional_table(frame)
    assert result.loc[0, "reaction_ec_numbers"] == ["2", "3"]
    assert bool(result.loc[0, "multifunctional_enzyme"]) is True


@pytest.mark.integration
def test_postprocess_target_table__matches_expected_output(
    tmp_path: Path,
    snapshot_resource: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "output.target_20250101.csv"
    input_csv = snapshot_resource / "target_organism_input.csv"
    source.write_bytes(input_csv.read_bytes())

    expected = pd.read_csv(
        snapshot_resource / "target_organism_expected.csv",
        dtype=object,
        keep_default_na=False,
    )

    Cellularity.get_lineage_names.cache_clear()

    def fake_lineage(cls, tax_id, email=None):
        if str(tax_id) == "2697049":
            return ("Viruses", "Betacoronavirus", "SARS-CoV-2")
        return ()

    monkeypatch.setattr(Cellularity, "get_lineage_names", classmethod(fake_lineage))

    output_path = Path(postprocess_target_table(str(source)))
    assert output_path.name == "organism.output.target_20250101.csv"

    result = pd.read_csv(output_path, dtype=object, keep_default_na=False)
    pdt.assert_frame_equal(result, expected)
