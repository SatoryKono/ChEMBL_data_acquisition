"""Unit tests for :mod:`library.postprocess.assays.steps`."""
from __future__ import annotations

import pandas as pd
import pytest

from library.postprocess.assays import steps


@pytest.mark.unit
def test_normalize_assay_metadata__respects_parameters() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["  primary  screen  "],
            "assay_test_type": ["Confirmatory"],
            "assay_format": [" cell based "],
        }
    )

    result = steps.normalize_assay_metadata(
        frame,
        uppercase_categories=False,
        strip_whitespace=False,
        collapse_internal_whitespace=True,
    )

    assert result.loc[0, "assay_type"] == " primary screen "
    assert result.loc[0, "assay_test_type"] == "Confirmatory"
    assert result.loc[0, "assay_format"] == " cell based "


@pytest.mark.unit
def test_enrich_assay_flags__custom_terms_and_default() -> None:
    frame = pd.DataFrame(
        {
            "assay_type": ["Primary screen", "Follow-up assay"],
        }
    )

    result = steps.enrich_assay_flags(
        frame,
        confirmatory_terms=("follow",),
        default_flag=True,
    )

    assert result["is_confirmatory"].tolist() == [True, True]


@pytest.mark.unit
def test_finalize_assay_records__enforce_schema_toggle(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, pd.DataFrame] = {}

    def _fake_validate(df: pd.DataFrame, *, context: str) -> pd.DataFrame:
        captured["frame"] = df.copy()
        captured["context"] = context
        return df.assign(validated=True)

    monkeypatch.setattr(steps, "validate_assays", _fake_validate)

    frame = pd.DataFrame(
        {
            "assay_chembl_id": [123],
            "assay_type": [" confirm "],
        }
    )

    enforced = steps.finalize_assay_records(frame, enforce_schema=True)
    assert captured["context"] == "assay_finalization"
    assert bool(enforced.loc[0, "validated"]) is True
    assert str(captured["frame"]["assay_chembl_id"].dtype) == "string"

    captured.clear()
    skipped = steps.finalize_assay_records(frame, enforce_schema=False, normalize_identifiers=False)
    assert "frame" not in captured
    assert skipped.equals(frame)
