from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from library.postprocessing.common.schema import DataFrameSchema, coerce_types


@pytest.mark.unit
def test_coerce_types_applies_schema_dtypes() -> None:
    frame = pd.DataFrame(
        {
            "symbol": ["BRCA1", "TP53"],
            "count": [1, None],
            "flag": [1, 0],
        }
    )

    schema = DataFrameSchema(
        dtypes={
            "symbol": "string",
            "count": "Int64",
            "flag": bool,
        }
    )

    result = coerce_types(frame, schema)

    expected = frame.copy()
    expected["symbol"] = expected["symbol"].astype("string")
    expected["count"] = expected["count"].astype("Int64")
    expected["flag"] = expected["flag"].astype(bool)

    pdt.assert_frame_equal(result, expected)

