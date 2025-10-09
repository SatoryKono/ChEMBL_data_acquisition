from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.target.postprocessing import postprocess_targets


@pytest.mark.unit
def test_postprocess_targets__uses_fallback_when_uniprot_id_missing() -> None:
    frame = pd.DataFrame(
        {
            "target_chembl_id": ["CHEMBL1"],
            "uniProtkbId": [pd.NA],
            "uniProtkbIdFallback": ["P12345-2"],
            "gene": ["GENE1|ALT1"],
            "geneName": [pd.NA],
        }
    )

    result = postprocess_targets(frame)

    assert result.loc[0, "uniprotkb_Id"] == "P12345"
