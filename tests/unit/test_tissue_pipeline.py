"""Unit tests for the tissue pipeline helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from library.pipelines.tissue import TISSUE_BASE_COLUMNS
from library.pipelines.tissue import pipeline as tissue_pipeline


@pytest.mark.parametrize(
    "limit, offset, expected",
    [
        (None, 0, ["CHEMBLT1", "CHEMBLT2", "#N/A", "CHEMBLT3"]),
        (2, 1, ["CHEMBLT2"]),
    ],
)
def test_read_tissue_identifiers__limit_and_offset(
    cfg,
    monkeypatch: pytest.MonkeyPatch,
    limit: int | None,
    offset: int,
    expected: list[str],
) -> None:
    """Ensure identifiers are trimmed and sliced deterministically."""

    calls: list[tuple[Path, str, object]] = []

    def fake_read_ids(path: Path, *, column: str, cfg: object):
        calls.append((path, column, cfg))
        yield from [
            " CHEMBLT1 ",
            " ",
            "CHEMBLT2",
            " #N/A ",
            "CHEMBLT3",
        ]

    monkeypatch.setattr(tissue_pipeline.io, "read_ids", fake_read_ids)

    identifiers = tissue_pipeline.read_tissue_identifiers(
        Path("dummy.csv"),
        column="tissue_chembl_id",
        io_cfg=cfg.io,
        limit=limit,
        offset=offset,
    )

    assert calls == [(Path("dummy.csv"), "tissue_chembl_id", cfg.io)]
    assert identifiers == expected


def test_prepare_tissue_dataframe__fills_missing_and_normalises() -> None:
    """Missing identifiers are appended and nullable fields normalised."""

    raw = pd.DataFrame(
        {
            "tissue_chembl_id": ["CHEMBLT1", "CHEMBLT1", "CHEMBLT2"],
            "pref_name": ["Alpha", "Alpha duplicate", ""],
            "uberon_id": ["UBERON:0001", "UBERON:0001", None],
            "irrelevant": ["ignored", "ignored", "ignored"],
        }
    )

    prepared, missing = tissue_pipeline.prepare_tissue_dataframe(
        raw,
        requested=["CHEMBLT1", "CHEMBLT3"],
    )

    assert missing == ("CHEMBLT3",)
    assert list(prepared.columns) == TISSUE_BASE_COLUMNS
    assert prepared["tissue_chembl_id"].tolist() == [
        "CHEMBLT1",
        "CHEMBLT2",
        "CHEMBLT3",
    ]
    assert (
        prepared.loc[prepared["tissue_chembl_id"] == "CHEMBLT2", "pref_name"]
        .isna()
        .all()
    )
    assert prepared.dtypes.tolist() == [pd.StringDtype()] * len(TISSUE_BASE_COLUMNS)
