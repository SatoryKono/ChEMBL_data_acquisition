from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pandas.testing as pdt
import pytest

from library.postprocessing.target.cellularity import add_cellularity_smart
from library.postprocessing.target.main import postprocess_target_table
from library.postprocessing.target.multifunctional import compute_multifunctional


@pytest.mark.integration
def test_postprocess_target_table__matches_power_query(tmp_path: Path, snapshot_resource: Path) -> None:
    source = snapshot_resource / "target_postprocess_input.csv"
    expected = snapshot_resource / "target_postprocess_expected.csv"
    working = tmp_path / source.name
    shutil.copy(source, working)

    output = Path(postprocess_target_table(working, fetcher=lambda tax_id, email: []))

    result = pd.read_csv(output, dtype=str, keep_default_na=False)
    expected_frame = pd.read_csv(expected, dtype=str, keep_default_na=False)
    pdt.assert_frame_equal(result, expected_frame)

    source_df = pd.read_csv(working, dtype=str, keep_default_na=False)
    base = source_df[
        [
            "target_chembl_id",
            "uniprot_id_primary",
            "organism",
            "taxon_id",
            "lineage_superkingdom",
            "lineage_phylum",
            "lineage_class",
        ]
    ].copy()
    for column in ("lineage_superkingdom", "lineage_phylum", "lineage_class"):
        base[column] = base[column].astype("string").str.lower()
    with_cellularity = add_cellularity_smart(
        base,
        "taxon_id",
        "lineage_superkingdom",
        "lineage_class",
        fetcher=lambda tax_id, email: [],
    )
    multifunctional = compute_multifunctional(source_df)
    joined = with_cellularity.merge(
        multifunctional[["target_chembl_id", "multifunctional_enzyme"]],
        on="target_chembl_id",
        how="left",
        sort=False,
    )
    assert joined["multifunctional_enzyme"].dtype == bool
