from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.target.postprocessing import finalise_targets
from library.schemas.targets import CELLULARITY_COLUMN_NAME


@pytest.mark.unit
def test_finalise_targets__derives_genus_from_organism_column() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniprotkb_Id": ["P12345"],
            "organism": ["Homo sapiens"],
            "lineage_superkingdom": ["Eukaryota"],
            "lineage_phylum": ["Chordata"],
            "lineage_class": ["Mammalia"],
        }
    )

    result = finalise_targets(frame)

    assert result.loc[0, "organism"] == "Homo"
    assert result.loc[0, CELLULARITY_COLUMN_NAME] == "Multicellular organism"
