from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from library.postprocessing import target as target_pp


def _make_classifier(**kwargs: object) -> target_pp.CellularityClassifier:
    return target_pp.CellularityClassifier(**kwargs)


def test_normalize_text__handles_none_and_trimming() -> None:
    assert target_pp.normalize_text(None) == ""
    assert target_pp.normalize_text("  AbC  ") == "abc"
    assert target_pp.normalize_text("β") == "β"


def test_first_element_text__returns_text_when_present() -> None:
    xml = textwrap.dedent(
        """
        <Taxon>
          <Lineage>cellular organisms</Lineage>
          <ScientificName>Example species</ScientificName>
        </Taxon>
        """
    ).strip()
    node = target_pp.ET.fromstring(xml)
    result = target_pp.first_element_text(list(node), "ScientificName")
    assert result == "Example species"


def test_first_element_text__returns_none_when_missing() -> None:
    xml = "<Taxon><Lineage>cellular organisms</Lineage></Taxon>"
    node = target_pp.ET.fromstring(xml)
    assert target_pp.first_element_text(list(node), "ScientificName") is None


@pytest.mark.parametrize(
    "superkingdom, phylum, expected",
    [
        ("Viruses", "", "acellular (virus)"),
        ("Bacteria", None, "unicellular"),
        ("Archaea", "", "unicellular"),
        ("Eukaryota", "Chordata", "multicellular"),
        ("Eukaryota", "ciliophora", "unicellular"),
        ("Eukaryota", "rhodophyta", "multicellular"),
        ("Eukaryota", "unknown", "ambiguous"),
    ],
)
def test_classify_by_lineage__scenarios(superkingdom: str, phylum: str | None, expected: str) -> None:
    classifier = _make_classifier()
    assert classifier.classify_by_lineage(superkingdom, phylum) == expected


def test_get_lineage_names__parses_lineage_and_scientific_name(monkeypatch: pytest.MonkeyPatch) -> None:
    xml = textwrap.dedent(
        """
        <TaxaSet>
          <Taxon>
            <Lineage>cellular organisms; Eukaryota; Metazoa; Chordata</Lineage>
            <ScientificName>Homo sapiens</ScientificName>
          </Taxon>
        </TaxaSet>
        """
    ).strip()

    def fake_http(url: str, params: dict[str, str], *, timeout: float) -> bytes:
        return xml.encode("utf-8")

    monkeypatch.setattr(target_pp, "_http_get", fake_http)
    classifier = _make_classifier()
    names = classifier._fetch_lineage_names("9606")
    assert names == [
        "cellular organisms",
        "Eukaryota",
        "Metazoa",
        "Chordata",
        "Homo sapiens",
    ]


def test_get_lineage_names__returns_empty_on_invalid_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_http(url: str, params: dict[str, str], *, timeout: float) -> bytes | None:
        return None

    monkeypatch.setattr(target_pp, "_http_get", fake_http)
    classifier = _make_classifier()
    assert classifier._fetch_lineage_names("0") == []


def test_classify_by_fetch__uses_lineage_names(monkeypatch: pytest.MonkeyPatch) -> None:
    classifier = _make_classifier()

    def fake_fetch(self: target_pp.CellularityClassifier, tax_id: str) -> list[str]:
        return ["Eukaryota", "Metazoa", "Chordata"]

    monkeypatch.setattr(target_pp.CellularityClassifier, "_fetch_lineage_names", fake_fetch)
    assert classifier.classify_by_fetch("9606") == "multicellular"


def test_add_cellularity_smart__lineage_without_fetch(monkeypatch: pytest.MonkeyPatch) -> None:
    classifier = _make_classifier()

    def forbid_fetch(_: object) -> str:
        raise AssertionError("fetch should not be called")

    monkeypatch.setattr(classifier, "classify_by_fetch", forbid_fetch)

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "taxon_id": [111],
            "lineage_superkingdom": ["bacteria"],
            "lineage_class": ["gammaproteobacteria"],
        }
    )
    result = classifier.add_cellularity_smart(
        frame,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
    )
    assert result.tolist() == ["unicellular"]


def test_add_cellularity_smart__fallback_fetch_when_ambiguous(monkeypatch: pytest.MonkeyPatch) -> None:
    classifier = _make_classifier()
    calls: list[str] = []

    def fake_fetch(tax_id: object) -> str:
        calls.append(str(tax_id))
        return "multicellular"

    monkeypatch.setattr(classifier, "classify_by_fetch", fake_fetch)

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL2"],
            "taxon_id": [9606],
            "lineage_superkingdom": [None],
            "lineage_class": [None],
        }
    )
    result = classifier.add_cellularity_smart(
        frame,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
    )
    assert result.tolist() == ["multicellular"]
    assert calls == ["9606"]


def test_multifunctional__reaction_ec_numbers_processing() -> None:
    values = target_pp._transform_reaction_ec_numbers("1.1.1.1|2.7.11.1|1.1.1.2|2.7.11.1")
    assert values == ["1", "2"]


def test_multifunctional__drops_columns_and_sets_flag() -> None:
    df = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL3"],
            "reaction_ec_numbers": ["1.1.1.1|2.7.11.1"],
            "isoform_ids": ["should be dropped"],
        }
    )
    result = target_pp._prepare_multifunctional_frame(df)
    assert list(result.columns) == ["target_chembl_id", "multifunctional_enzyme"]
    assert bool(result.loc[0, "multifunctional_enzyme"])


def test_process_targets__produces_expected_columns_and_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    data = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1", "CHEMBL2"],
            "uniprot_id_primary": ["P12345", "Q54321"],
            "organism": ["Species A", "Species B"],
            "taxon_id": [111, 222],
            "lineage_superkingdom": ["viruses", None],
            "lineage_phylum": ["", ""],
            "lineage_class": ["", ""],
            "reaction_ec_numbers": ["1.1.1.1|2.7.11.1", "1.1.1.1"],
        }
    )
    input_path = tmp_path / "output.target_20251004.csv"
    data.to_csv(input_path, index=False)

    def forbid_http(url: str, params: dict[str, str], *, timeout: float) -> bytes:
        raise AssertionError("HTTP should not be called in offline mode")

    monkeypatch.setattr(target_pp, "_http_get", forbid_http)

    output_path = target_pp.process_targets(
        input_csv=str(input_path),
        output_prefix="organism",
        offline=True,
        verbose=False,
    )

    assert output_path.name == "organism.output.target_20251004.csv"
    result = pd.read_csv(output_path)
    assert list(result.columns) == list(target_pp.TARGET_OUTPUT_COLUMNS)
    assert result.loc[0, "cellularity"] == "acellular (virus)"
    assert result.loc[1, "cellularity"] == "ambiguous"
    assert bool(result.loc[0, "multifunctional_enzyme"])
    assert not bool(result.loc[1, "multifunctional_enzyme"])


def test_process_targets__offline_ambiguous_without_http(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL10"],
            "uniprot_id_primary": ["P00010"],
            "organism": ["Organism"],
            "taxon_id": [333],
            "lineage_superkingdom": [None],
            "lineage_phylum": [None],
            "lineage_class": [None],
            "reaction_ec_numbers": ["1.1.1.1|3.2.1.4"],
        }
    )
    input_path = tmp_path / "output.target_20250101.csv"
    frame.to_csv(input_path, index=False)

    def forbid_http(url: str, params: dict[str, str], *, timeout: float) -> bytes:
        raise AssertionError("HTTP should not be executed in offline mode")

    monkeypatch.setattr(target_pp, "_http_get", forbid_http)

    output_path = target_pp.process_targets(
        input_csv=str(input_path),
        output_prefix="organism",
        offline=True,
        verbose=False,
    )
    df = pd.read_csv(output_path)
    assert df.loc[0, "cellularity"] == "ambiguous"


def test_process_targets__fetch_classification(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    xml = textwrap.dedent(
        """
        <TaxaSet>
          <Taxon>
            <Lineage>cellular organisms; Eukaryota; Metazoa; Chordata</Lineage>
            <ScientificName>Example organism</ScientificName>
          </Taxon>
        </TaxaSet>
        """
    ).strip()

    def fake_http(url: str, params: dict[str, str], *, timeout: float) -> bytes:
        return xml.encode("utf-8")

    monkeypatch.setattr(target_pp, "_http_get", fake_http)

    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL20"],
            "uniprot_id_primary": ["P00020"],
            "organism": ["Organism"],
            "taxon_id": [9606],
            "lineage_superkingdom": ["Eukaryota"],
            "lineage_phylum": ["Chordata"],
            "lineage_class": ["Mammalia"],
            "reaction_ec_numbers": ["1.1.1.1|2.7.11.1"],
        }
    )
    input_path = tmp_path / "output.target_20250102.csv"
    frame.to_csv(input_path, index=False)

    output_path = target_pp.process_targets(
        input_csv=str(input_path),
        output_prefix="organism",
        offline=False,
        verbose=False,
    )

    df = pd.read_csv(output_path)
    assert df.loc[0, "cellularity"] == "multicellular"

