from __future__ import annotations

from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.config import Config
from library.pipelines.target import cellularity, helpers, multifunctional
from library.postprocessing.target.cellularity import Cellularity
from library.schemas.targets import CELLULARITY_COLUMN_NAME, TARGETS_COLUMN_ORDER

INPUT_FILE = "target_postprocess_power_query_input.csv"
EXPECTED_FILE = "target_postprocess_power_query_expected.csv"


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
    assert working_output.read_bytes() == expected_path.read_bytes()

    result = pd.read_csv(working_output, dtype=str, keep_default_na=False).astype(
        "string"
    )
    expected = pd.read_csv(expected_path, dtype=str, keep_default_na=False).astype(
        "string"
    )
    pdt.assert_frame_equal(result, expected)


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
