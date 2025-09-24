"""Behavioural checks for :func:`library.uniprot_library.process`."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from pytest import MonkeyPatch

from library import uniprot_library as ul
from library.config import UniprotCfg


def test_process_writes_expected_columns(
    tmp_path: Path, monkeypatch: MonkeyPatch
) -> None:
    """The UniProt CSV writer respects :data:`UNIPROT_OUTPUT_COLUMNS`."""

    input_csv = tmp_path / "ids.csv"
    input_csv.write_text("uniprot_id\nP12345\n", encoding="utf-8")
    output_csv = tmp_path / "out.csv"

    sample = {column: f"{column}_value" for column in ul.UNIPROT_OUTPUT_COLUMNS}
    sample["uniprot_id"] = "P12345"

    def fake_collect_info(
        uid: str, data_dir: Path | str | None = None, *, cfg: UniprotCfg
    ) -> dict[str, object]:
        assert uid == "P12345"
        assert isinstance(cfg, UniprotCfg)
        payload = {**sample, "secondaryAccessions": ["S1", "S2"], "extra": "x"}
        return payload

    monkeypatch.setattr(ul, "collect_info", fake_collect_info)

    ul.process(
        str(input_csv),
        str(output_csv),
        data_dir=tmp_path / "json",
        cfg=UniprotCfg(),
    )

    df = pd.read_csv(output_csv, dtype=str)
    assert df.columns.tolist() == ul.UNIPROT_OUTPUT_COLUMNS
    row = df.loc[0]
    assert row["uniprot_id"] == "P12345"
    assert row["secondaryAccessions"] == "S1|S2"

