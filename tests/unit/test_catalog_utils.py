from __future__ import annotations

import pandas as pd
import pytest

from library.pipelines.testitem import catalog


@pytest.mark.unit
def test_ensure_no_parant_column__raises() -> None:
    frame = pd.DataFrame({"parant_molecule_id": ["CHEMBL1"]})
    with pytest.raises(ValueError):
        catalog.ensure_no_parant_column(frame)


@pytest.mark.unit
@pytest.mark.parametrize(
    "value, child_id, expected",
    [
        (None, "CHEMBL1", None),
        ("", "CHEMBL1", None),
        ("NULL", "CHEMBL1", None),
        ("CHEMBL1", "CHEMBL1", None),
        ("NO PARENT", "CHEMBL2", None),
        (" chembl2 ", "CHEMBL1", "CHEMBL2"),
    ],
)
def test_normalise_parent_identifier__cases(
    value: object, child_id: str, expected: str | None
) -> None:
    assert catalog._normalise_parent_identifier(value, child_id=child_id) == expected


@pytest.mark.unit
def test_merge_parent_stats__accumulates_counts() -> None:
    base = catalog.ParentLookupStats(
        source=catalog.PARENT_LOOKUP_SOURCE_CACHE,
        missing=1,
        unique=2,
        attached=1,
        uncovered=0,
        failed_ids=("CHEMBL1",),
        hierarchy_attached=1,
        fallback_attached=0,
        no_parent=0,
    )
    update = catalog.ParentLookupStats(
        source=catalog.PARENT_LOOKUP_SOURCE_LOOKUP,
        missing=2,
        unique=1,
        attached=1,
        uncovered=1,
        failed_ids=("CHEMBL2", "CHEMBL1"),
        hierarchy_attached=0,
        fallback_attached=1,
        no_parent=1,
    )

    merged = catalog._merge_parent_stats(base, update)

    assert merged.source == catalog.PARENT_LOOKUP_SOURCE_LOOKUP
    assert merged.missing == 3
    assert merged.unique == 3
    assert merged.attached == 2
    assert merged.uncovered == 1
    assert merged.hierarchy_attached == 1
    assert merged.fallback_attached == 1
    assert merged.no_parent == 1
    assert merged.failed_ids == ("CHEMBL1", "CHEMBL2")


@pytest.mark.unit
def test_merge_parent_stats__keeps_higher_priority_source() -> None:
    base = catalog.ParentLookupStats(
        source=catalog.PARENT_LOOKUP_SOURCE_SYNC,
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )
    update = catalog.ParentLookupStats(
        source=catalog.PARENT_LOOKUP_SOURCE_CACHE,
        missing=0,
        unique=0,
        attached=0,
        uncovered=0,
    )

    merged = catalog._merge_parent_stats(base, update)
    assert merged.source == catalog.PARENT_LOOKUP_SOURCE_SYNC
