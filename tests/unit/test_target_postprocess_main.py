from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt

import pytest

from library.postprocessing.target.main import postprocess_target_table

@pytest.mark.unit  # type: ignore[misc]
def test_postprocess_target_table__produces_power_query_equivalent_csv(
    tmp_path: Path, snapshot_resource: Path
) -> None:
    input_path = snapshot_resource / "target_postprocess_input.csv"
    expected_path = snapshot_resource / "target_postprocess_expected.csv"

    working_input = tmp_path / input_path.name
    working_input.write_bytes(input_path.read_bytes())

    output_location = postprocess_target_table(working_input)
    output_path = Path(output_location)

    assert output_path == tmp_path / f"organism.{input_path.name}"
    # Normalise EOL for cross-platform determinism
    assert output_path.read_bytes().replace(b"\r\n", b"\n") == expected_path.read_bytes().replace(b"\r\n", b"\n")

    result_frame = pd.read_csv(output_path, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected_path, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result_frame, expected_frame)

    bool_frame = pd.read_csv(
        output_path,
        dtype={"multifunctional_enzyme": "boolean"},
        keep_default_na=False,
    )
    assert bool_frame["multifunctional_enzyme"].dtype.name == "boolean"


@pytest.mark.unit  # type: ignore[misc]
def test_postprocess_target_table__fills_missing_columns(tmp_path: Path) -> None:
    input_path = tmp_path / "targets_minimal.csv"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "taxon_id": ["9606"],
            "reaction_ec_numbers": ["1.1.1.1|2.2.2.2"],
        }
    )
    frame.to_csv(input_path, index=False)

    fetch_calls: list[tuple[object, object | None]] = []

    def _fake_fetcher(tax_id: object, email: str | None) -> list[str]:
        fetch_calls.append((tax_id, email))
        return ["Eukaryota", "Chordata"]

    output_location = postprocess_target_table(input_path, fetcher=_fake_fetcher)
    output_path = Path(output_location)

    result_frame = pd.read_csv(output_path, dtype=str, keep_default_na=False)

    for column in [
        "uniprot_id_primary",
        "organism",
        "lineage_superkingdom",
        "lineage_phylum",
        "lineage_class",
    ]:
        assert column in result_frame.columns
        assert result_frame.at[0, column] == ""

    assert fetch_calls == [("9606", None)]
    assert result_frame.at[0, "cellularity"] == "multicellular"
    assert result_frame.at[0, "multifunctional_enzyme"] == "True"
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprot_id_primary": ["P11111"],
            "organism": ["Human"],
            "taxon_id": ["9606"],
            "lineage_superkingdom": ["Eukaryota"],
            "lineage_phylum": ["Chordata"],
            "lineage_class": ["Mammalia"],
        }
    )
    frame.to_csv(input_path, index=False)

    output_location = postprocess_target_table(input_path)
    output_path = Path(output_location)

    assert output_path.name == "organism.output.target_20250101.csv"
    assert output_path.exists()

@pytest.mark.unit  # type: ignore[misc]
def test_postprocess_target_table__handles_missing_identifier_column(
    tmp_path: Path,
) -> None:
    input_path = tmp_path / "targets_missing_id.csv"
    frame = pd.DataFrame(
        {
            "taxon_id": ["9606"],
            "reaction_ec_numbers": ["1.2.3.4|2.7.11.1"],
        }
    )
    frame.to_csv(input_path, index=False)

    fetch_calls: list[tuple[object, object | None]] = []

    def _fake_fetcher(tax_id: object, email: str | None) -> list[str]:
        fetch_calls.append((tax_id, email))
        return ["Eukaryota", "Chordata"]

    output_location = postprocess_target_table(input_path, fetcher=_fake_fetcher)
    output_path = Path(output_location)

    result_frame = pd.read_csv(output_path, dtype=str, keep_default_na=False)

    assert result_frame.at[0, "target_chembl_id"] == ""
    assert result_frame.at[0, "multifunctional_enzyme"] == "True"
    assert fetch_calls == [("9606", None)]

